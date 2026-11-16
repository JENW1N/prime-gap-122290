"""Test a presieved search with CGBN, checkpointing completed endpoint pairs."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
from time import perf_counter

import gmpy2

from gpu_fermat import FermatGpu


def store(path, data):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="ascii")
    temporary.replace(path)


def load_search(folder):
    meta = json.loads((folder / "parameters.json").read_text())
    rows = []
    expected_header = f"ABC {meta['expr']} & {meta['expr']}+{meta['gap']}"
    for shard in sorted(folder.glob("shard*")):
        with (shard / "candidates.txt").open(encoding="ascii") as stream:
            if stream.readline().rstrip() != expected_header:
                raise ValueError(f"Inconsistent expression in {shard.name}")
            rows.extend(map(int, stream))
    first, count = meta["first"], meta["rows"]
    if (len(rows) != meta["eligible"] or any(a >= b for a, b in zip(rows, rows[1:]))
            or not first <= rows[0] <= rows[-1] < first + count):
        raise ValueError("Invalid presieved candidate list")
    return meta, rows


def reference_passes(shard):
    hits = []
    for line in (shard / "stdout.log").read_bytes().splitlines():
        if b"is 2-PRP!" in line:
            match = re.match(rb"(\d+)\*", line.strip())
            if match is None:
                raise ValueError("Unrecognized PFGW result")
            hits.append(int(match[1]))
    if len(hits) != len(set(hits)):
        raise ValueError("Reference shard contains a double probable-prime row")
    return set(hits)


def independently_check(values, residues, flags):
    for index, n in enumerate(values):
        exact = gmpy2.powmod(2, n - 1, n)
        if residues[index] != int(exact % (1 << 64)) or flags[index] != (exact == 1):
            raise ArithmeticError(f"GPU residue check failed at item {index}")


def prove_partial_interior(left, gap, uncovered, gpu=None):
    from search_exact_gap import get_primes

    remaining = set(uncovered)
    witnesses = []
    for prime in get_primes(200_000):
        q = int(prime)
        if q >= left:
            break
        first = int((-left) % q)
        for offset in range(first, gap, q):
            if offset in remaining:
                remaining.remove(offset)
                witnesses.append([offset, "divisor", q])

    offsets = sorted(remaining)
    if gpu is not None and offsets:
        _, flags, _ = gpu.run([left + offset for offset in offsets])
        for offset, passed in zip(offsets, flags):
            if passed and gmpy2.is_probab_prime(left + offset, 12):
                print(f"Rejected the pair because offset {offset} is probably prime.", flush=True)
                return None

    for offset in offsets:
        number = left + offset
        for base in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
            if gmpy2.powmod(base, number - 1, number) != 1:
                witnesses.append([offset, "fermat", base])
                break
            if base == 2 and gmpy2.is_probab_prime(number, 12):
                print(f"Rejected the pair because offset {offset} is probably prime.", flush=True)
                return None
        else:
            print(f"Rejected the pair because offset {offset} lacks a compositeness witness.", flush=True)
            return None
    return sorted(witnesses)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("--batch", type=int, default=1024)
    parser.add_argument("--checkpoint-batches", type=int, default=8)
    parser.add_argument("--reference-shard", type=str)
    args = parser.parse_args()
    meta, rows = load_search(args.folder)
    a, modulus = gmpy2.mpz(meta["a"]), gmpy2.mpz(meta["step"])
    gap = meta["gap"]
    if args.reference_shard:
        shard = args.folder / args.reference_shard
        rows = list(map(int, (shard / "candidates.txt").read_text().splitlines()[1:]))
        expected = reference_passes(shard)
    bits = ((a + rows[-1] * modulus + gap).bit_length() // 1024 + 1) * 1024
    if bits > 32768:
        raise ValueError("Candidates exceed the CGBN width limit")
    gpu = FermatGpu(bits, 32, 512)
    uncovered = meta.get("uncovered", [])
    identity_data = (meta["a"] + "\n" + meta["step"] + "\n" + str(gap) + "\n"
                     + "".join(f"{row}\n" for row in rows))
    if uncovered:
        identity_data += "interior offsets\n" + "".join(f"{offset}\n" for offset in uncovered)
    identity = sha256(identity_data.encode("ascii")).hexdigest()
    state_path = args.folder / "gpu_checkpoint.json"
    candidate_path = args.folder / "gpu_candidate.json"
    state = {"input_sha256": identity, "tested": 0, "left_prps": 0, "double_base2": []}
    if not args.reference_shard and state_path.exists():
        state = json.loads(state_path.read_text())
        if state["input_sha256"] != identity:
            raise ValueError("Checkpoint belongs to different rows")
    if not args.reference_shard and candidate_path.exists():
        raise ValueError(f"Inspect existing candidate in {candidate_path.name}")
    starting = state["tested"]
    if not 0 <= starting <= len(rows):
        raise ValueError("Checkpoint outside the input")
    print(f"{len(rows)} presieved rows, starting at offset {starting}, {bits}-bit kernel", flush=True)

    left_passes = []
    begin = perf_counter()
    for index in range(starting, len(rows), args.batch):
        batch = rows[index:index + args.batch]
        values = [a + m * modulus for m in batch]
        residues, flags, _ = gpu.run(values)
        if index == starting:
            independently_check(values[:2], residues[:2], flags[:2])
        left_passes.extend(row for row, passed in zip(batch, flags) if passed)
        tested = index + len(batch)
        if args.reference_shard:
            continue
        batch_number = (index - starting) // args.batch + 1
        if batch_number % args.checkpoint_batches and tested != len(rows):
            continue

        state["tested"] = tested
        state["left_prps"] += len(left_passes)
        if left_passes:
            right_values = [a + row * modulus + gap for row in left_passes]
            _, flags, _ = gpu.run(right_values)
            for row, right, passed in zip(left_passes, right_values, flags):
                if passed:
                    state["double_base2"].append(row)
                    left = right - gap
                    if gmpy2.is_probab_prime(left, 12) and gmpy2.is_probab_prime(right, 12):
                        candidate = {
                            "row": row, "gap": gap, "left": str(left), "right": str(right),
                            "input_sha256": identity
                        }
                        if uncovered:
                            print(f"Testing the partial cover at double probable-prime row {row}.", flush=True)
                            interior = prove_partial_interior(left, gap, uncovered, gpu)
                            if interior is None:
                                continue
                            candidate["interior_composites"] = interior
                        store(candidate_path, candidate)
                        store(state_path, state)
                        print(f"Found a double probable-prime candidate at row {row}.", flush=True)
                        return
            left_passes.clear()
        store(state_path, state)
        seconds = perf_counter() - begin
        print(f"Tested {tested}/{len(rows)} rows, {state['left_prps']} left passes, "
              f"{(tested-starting)/seconds:.1f} rows/s", flush=True)

    if args.reference_shard:
        if set(left_passes) != expected:
            raise ArithmeticError(f"GPU/PFGW mismatch at rows {sorted(set(left_passes) ^ expected)[:20]}")
        print(f"Agreed with PFGW on all {len(rows)} rows, including {len(expected)} probable primes, "
              f"in {perf_counter()-begin:.2f}s.")
    else:
        print("Completed all rows without a double probable-prime candidate.")


if __name__ == "__main__":
    main()
