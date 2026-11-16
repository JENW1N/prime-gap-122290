"""Build and search a CRT progression for a specified consecutive-prime gap."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys
from time import perf_counter

import gmpy2
import numpy as np
from numba import njit


sys.set_int_max_str_digits(0)


@njit
def sieve_rows(size, first, primes, residues1, residues2):
    allowed = np.ones(size, dtype=np.bool_)
    for i in range(len(primes)):
        p = primes[i]
        for j in range((residues1[i] - first) % p, size, p):
            allowed[j] = False
        for j in range((residues2[i] - first) % p, size, p):
            allowed[j] = False
    return allowed


def get_primes(limit):
    sieve = np.ones(limit + 1, dtype=np.bool_)
    sieve[:2] = False
    for q in range(2, math.isqrt(limit) + 1):
        if sieve[q]:
            sieve[q * q::q] = False
    return np.flatnonzero(sieve)


def make_progression(mask):
    h = mask["gap"]
    residues = {2: 1}
    residues.update({int(q): (-2 * r) % int(q) for q, r in mask["offsets"].items()})
    uncovered = np.ones(h + 1, dtype=np.bool_)
    a, step = gmpy2.mpz(0), gmpy2.mpz(1)
    for q, b in sorted(residues.items()):
        if b == 0 or (b + h) % q == 0:
            raise ArithmeticError(f"Endpoint hit by prime {q}")
        uncovered[(-b) % q::q] = False
        a += step * (((b - a) * gmpy2.invert(step, q)) % q)
        step *= q
    actual = np.flatnonzero(uncovered).tolist()
    if actual != [0, *mask["uncovered"], h]:
        raise ArithmeticError("Mask survivors do not match its assignments")
    return a, step


def prepare(mask_path, target, first, rows, sieve_bound, shards):
    begin = perf_counter()
    mask = json.loads(mask_path.read_text())
    h = mask["gap"]
    a, step = make_progression(mask)
    mod_primes = set(map(int, mask["offsets"])) | {2}
    primes = np.array([int(q) for q in get_primes(sieve_bound) if q not in mod_primes], dtype=np.int64)
    residues1 = np.empty(len(primes), dtype=np.int64)
    residues2 = np.empty(len(primes), dtype=np.int64)
    for i, q64 in enumerate(primes):
        q = int(q64)
        inverse = int(gmpy2.invert(step, q))
        left = int(a % q)
        residues1[i] = (-left * inverse) % q
        residues2[i] = ((-left - h) * inverse) % q
    allowed = sieve_rows(rows, first, primes, residues1, residues2)
    eligible = np.flatnonzero(allowed) + first
    largest = max(mod_primes)
    quotient = math.prod(int(q) for q in get_primes(largest) if q not in mod_primes)
    expr = f"$a*({largest}#/{quotient})+{a}"
    target.mkdir(parents=True, exist_ok=True)
    run = {
        "gap": h,
        "a": str(a),
        "step": str(step),
        "uncovered": mask["uncovered"],
        "first": first,
        "rows": rows,
        "sieve_bound": sieve_bound,
        "eligible": len(eligible),
        "expr": expr,
    }
    for index, values in enumerate(np.array_split(eligible, shards)):
        directory = target / f"shard{index:02d}"
        directory.mkdir(exist_ok=True)
        (directory / "candidates.txt").write_text(
            f"ABC {expr} & {expr}+{h}\n" + "".join(f"{row}\n" for row in values),
            encoding="ascii",
        )
    (target / "parameters.json").write_text(json.dumps(run))
    print(
        f"{h=}, {len(str(step))} digits, {len(mask['uncovered'])} interior survivors, "
        f"{len(eligible)} of {rows} rows passed presieve to {sieve_bound} "
        f"in {perf_counter() - begin:.2f}s",
        flush=True,
    )


def run_shards(target, engine):
    processes = []
    start = perf_counter()
    try:
        for folder in sorted(target.glob("shard*")):
            log = (folder / "stdout.log").open("w", encoding="ascii", errors="replace")
            child = subprocess.Popen(
                [str(engine), "candidates.txt", "-b2", "-f0", "-u0", "-Cnormal", "-N"],
                cwd=folder,
                stdout=log,
                stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW,
            )
            processes.append((folder, child, log))
        failed = []
        for folder, child, log in processes:
            code = child.wait()
            log.close()
            if code != 0:
                failed.append((folder.name, code))
        if failed:
            raise RuntimeError(f"PFGW failures {failed}")
    finally:
        for _, child, log in processes:
            if child.poll() is None:
                child.terminate()
                child.wait()
            log.close()
    print(f"Completed {len(processes)} PFGW shards in {perf_counter() - start:.2f}s", flush=True)


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="action", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("mask", type=Path)
    prep.add_argument("target", type=Path)
    prep.add_argument("--first", type=int, default=1)
    prep.add_argument("--rows", type=int, default=100000)
    prep.add_argument("--sieve", type=int, default=10**7)
    prep.add_argument("--shards", type=int, default=8)
    search = sub.add_parser("run")
    search.add_argument("target", type=Path)
    search.add_argument("--engine", type=Path, default=Path("work/pfgw/distribution/pfgw64.exe"))
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args.mask, args.target, args.first, args.rows, args.sieve, args.shards)
    else:
        run_shards(args.target, args.engine.resolve(strict=True))


if __name__ == "__main__":
    main()
