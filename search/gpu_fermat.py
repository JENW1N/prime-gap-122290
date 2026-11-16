"""Benchmark and run chunked CGBN base-two Fermat tests on CRT candidates."""

import argparse
import ctypes
import hashlib
import json
import os
from pathlib import Path
import sys
from time import perf_counter

import gmpy2
import numpy as np


WORK = Path(__file__).resolve().parent
DEPS = WORK / "gpu_deps"
sys.path.insert(0, str(DEPS))
NVRTC_DIRECTORY = DEPS / "nvidia" / "cuda_nvrtc" / "bin"
CUDA_DLLS = os.add_dll_directory(str(NVRTC_DIRECTORY))
CUDA_COMPILER = ctypes.CDLL(str(NVRTC_DIRECTORY / "nvrtc64_120_0.dll"))
from cuda.bindings import driver, nvrtc


def result(call):
    status, *values = call
    if int(status):
        raise RuntimeError(str(status))
    if len(values) == 1:
        return values[0]
    return tuple(values)


KERNEL = r"""
typedef unsigned int uint32_t;
typedef int int32_t;
typedef unsigned long long uint64_t;
typedef long long int64_t;
extern "C" __device__ int printf(const char *, ...);
#include "cgbn_nvrtc.h"

typedef cgbn_context_t<GROUP_SIZE> context_t;
typedef cgbn_env_t<context_t, WORD_BITS> env_t;
typedef cgbn_mem_t<WORD_BITS> words_t;

extern "C" __global__ void iterate(words_t *moduli, words_t *states,
                                   const uint32_t *np0, int count, int hi, int lo) {
    int index = (blockDim.x * blockIdx.x + threadIdx.x) / GROUP_SIZE;
    if (index >= count) return;
    context_t context(cgbn_no_checks, nullptr, index);
    env_t env(context);
    env_t::cgbn_t n, r;
    cgbn_load(env, n, moduli + index);
    cgbn_load(env, r, states + index);
    uint32_t inv = np0[index];
    #pragma unroll 1
    for (int bit = hi; bit >= lo; --bit) {
        cgbn_mont_sqr(env, r, r, n, inv);
        if (bit > 0 && cgbn_extract_bits_ui32(env, n, bit, 1)) {
            cgbn_add(env, r, r, r);
            if (cgbn_compare(env, r, n) >= 0) cgbn_sub(env, r, r, n);
        }
    }
    cgbn_store(env, states + index, r);
}

extern "C" __global__ void finish(words_t *moduli, words_t *states, const uint32_t *np0,
                                  uint64_t *res64, uint32_t *flags, int count) {
    int index = (blockDim.x * blockIdx.x + threadIdx.x) / GROUP_SIZE;
    if (index >= count) return;
    context_t context(cgbn_no_checks, nullptr, index);
    env_t env(context);
    env_t::cgbn_t n, r;
    cgbn_load(env, n, moduli + index);
    cgbn_load(env, r, states + index);
    cgbn_mont2bn(env, r, r, n, np0[index]);
    uint32_t low = cgbn_extract_bits_ui32(env, r, 0, 32);
    uint32_t high = cgbn_extract_bits_ui32(env, r, 32, 32);
    uint32_t pass = cgbn_compare_ui32(env, r, 1) == 0;
    if (threadIdx.x % GROUP_SIZE == 0) {
        res64[index] = low | ((uint64_t) high << 32);
        flags[index] = pass;
    }
}
"""


def compile_kernel(bits, tpi, kernel=KERNEL):
    include = WORK / "CGBN" / "include" / "cgbn"
    header = (include / "cgbn.h").read_text(encoding="ascii")
    header = "\n".join(line for line in header.splitlines() if not line.startswith(
        ('#include <', '#include "cgbn.cu"', 'cudaError_t cgbn_', 'bool        cgbn_error_',
         'void        cgbn_error_', 'const char *cgbn_error_')
    ))
    code = kernel.encode("ascii")
    stamp = hashlib.sha256(code + header.encode("ascii") + str((bits, tpi)).encode()).hexdigest()[:16]
    cache = WORK / f"fermat_{stamp}.cubin"
    if cache.exists():
        return cache.read_bytes()
    options = [b"--gpu-architecture=sm_120", b"--std=c++17", b"--device-as-default-execution-space",
               f"-I{include}".encode(), f"-DGROUP_SIZE={tpi}".encode(), f"-DWORD_BITS={bits}".encode()]
    program = result(nvrtc.nvrtcCreateProgram(code, b"fermat.cu", 1, [header.encode("ascii")], [b"cgbn_nvrtc.h"]))
    try:
        error, = nvrtc.nvrtcCompileProgram(program, len(options), options)
        logsize = result(nvrtc.nvrtcGetProgramLogSize(program))
        log = bytearray(logsize)
        result(nvrtc.nvrtcGetProgramLog(program, log))
        if int(error):
            raise RuntimeError(log.decode(errors="replace").rstrip("\0"))
        size = result(nvrtc.nvrtcGetCUBINSize(program))
        cubin = bytearray(size)
        result(nvrtc.nvrtcGetCUBIN(program, cubin))
    finally:
        result(nvrtc.nvrtcDestroyProgram(program))
    cache.write_bytes(cubin)
    return bytes(cubin)


class FermatGpu:
    def __init__(self, bits, tpi, chunk_bits):
        self.bits, self.tpi, self.chunk_bits = bits, tpi, chunk_bits
        result(driver.cuInit(0))
        device = result(driver.cuDeviceGet(0))
        self.device = device
        self.context = result(driver.cuDevicePrimaryCtxRetain(device))
        result(driver.cuCtxSetCurrent(self.context))
        self.module = result(driver.cuModuleLoadData(compile_kernel(bits, tpi)))
        self.iterate = result(driver.cuModuleGetFunction(self.module, b"iterate"))
        self.finish = result(driver.cuModuleGetFunction(self.module, b"finish"))

    def run(self, values):
        count = len(values)
        assert count and all(n & 1 and n.bit_length() < self.bits for n in values)
        width = self.bits // 8
        mont_r = gmpy2.mpz(1) << self.bits
        moduli = np.frombuffer(b"".join(n.to_bytes(width, "little") for n in values), dtype=np.uint32)
        states = np.frombuffer(b"".join((mont_r % n).to_bytes(width, "little") for n in values), dtype=np.uint32)
        np0 = np.array([(-pow(int(n % (1 << 32)), -1, 1 << 32)) % (1 << 32) for n in values], dtype=np.uint32)
        res64 = np.empty(count, dtype=np.uint64)
        flags = np.empty(count, dtype=np.uint32)
        arrays = (moduli, states, np0, res64, flags)
        buffers = []
        begin = perf_counter()
        try:
            for i, array in enumerate(arrays):
                buffer = result(driver.cuMemAlloc(array.nbytes))
                buffers.append(buffer)
                if i < 3:
                    result(driver.cuMemcpyHtoD(buffer, array.ctypes.data, array.nbytes))
            threads = 128
            blocks = (count * self.tpi + threads - 1) // threads
            ptr, integer = ctypes.c_uint64, ctypes.c_int
            hi = max(n.bit_length() for n in values) - 1
            while hi >= 0:
                lo = max(0, hi - self.chunk_bits + 1)
                args = (tuple(map(int, buffers[:3])) + (count, hi, lo), (ptr, ptr, ptr, integer, integer, integer))
                result(driver.cuLaunchKernel(self.iterate, blocks, 1, 1, threads, 1, 1, 0, 0, args, 0))
                result(driver.cuCtxSynchronize())
                hi = lo - 1
            args = (tuple(map(int, buffers)) + (count,), (ptr, ptr, ptr, ptr, ptr, integer))
            result(driver.cuLaunchKernel(self.finish, blocks, 1, 1, threads, 1, 1, 0, 0, args, 0))
            for array, buffer in zip(arrays[3:], buffers[3:]):
                result(driver.cuMemcpyDtoH(array.ctypes.data, buffer, array.nbytes))
        finally:
            for buffer in buffers:
                result(driver.cuMemFree(buffer))
        return res64, flags, perf_counter() - begin


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search", type=Path)
    parser.add_argument("--count", type=int, default=8)
    parser.add_argument("--tpi", type=int, default=32)
    parser.add_argument("--chunk-bits", type=int, default=512)
    parser.add_argument("--verify", type=int, default=3)
    args = parser.parse_args()
    if args.search:
        meta = json.loads((args.search / "parameters.json").read_text())
        a, W = gmpy2.mpz(meta["a"]), gmpy2.mpz(meta["step"])
        first_shard = min(args.search.glob("shard*"))
        rows = list(map(int, (first_shard / "candidates.txt").read_text().splitlines()[1:1+args.count]))
        values = [a + m * W for m in rows]
    else:
        values = [gmpy2.mpz(n) for n in (561, 104729, 4294967311, 12345678901234567891)]
    group_bits = 32 * args.tpi
    bits = (max(n.bit_length() for n in values) // group_bits + 1) * group_bits
    print(f"Preparing {bits}-bit Fermat kernel on {args.tpi} threads per candidate", flush=True)
    gpu = FermatGpu(bits, args.tpi, args.chunk_bits)
    residues, flags, seconds = gpu.run(values)
    print(f"{len(values)} base-two tests in {seconds:.3f}s, {len(values)/seconds:.1f}/s", flush=True)
    for i, n in enumerate(values[:args.verify]):
        exact = gmpy2.powmod(2, n - 1, n)
        if residues[i] != int(exact % (1 << 64)) or flags[i] != (exact == 1):
            raise ArithmeticError(f"Incorrect GPU result at item {i}: {residues[i]} vs {int(exact % (1 << 64))}")
    print(f"Verified {min(args.verify, len(values))} against GMP. Base-two passes: {int(flags.sum())}.")


if __name__ == "__main__":
    main()
