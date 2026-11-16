"""Search CRT assignments while keeping the two endpoints available."""

import argparse
import json
import math
from pathlib import Path
from time import perf_counter

import numpy as np
from numba import njit
from sympy import primerange


@njit
def change(covers, uncovered, index, count, prime, residue, delta):
    if residue < 0:
        return count
    for site in range(residue, len(covers), prime):
        if delta < 0:
            covers[site] -= 1
            if covers[site] == 0:
                uncovered[count] = site
                index[site] = count
                count += 1
        else:
            if covers[site] == 0:
                count -= 1
                moved = uncovered[count]
                uncovered[index[site]] = moved
                index[moved] = index[site]
            covers[site] += 1
    return count


@njit
def best_residue(uncovered, count, prime, right, scores, touched):
    used = 0
    right_residue = right % prime
    for idx in range(count):
        residue = uncovered[idx] % prime
        if residue != 0 and residue != right_residue:
            if scores[residue] == 0:
                touched[used] = residue
                used += 1
            scores[residue] += 1
    best = -1
    value = 0
    ties = 0
    for idx in range(used):
        residue = touched[idx]
        score = scores[residue]
        scores[residue] = 0
        if score > value:
            best, value, ties = residue, score, 1
        elif score == value:
            ties += 1
            if np.random.randint(ties) == 0:
                best = residue
    return best


@njit
def optimize(primes, right, start, rounds, sweeps, seed):
    np.random.seed(seed)
    residues = start.copy()
    covers = np.zeros(right, dtype=np.int16)
    covers[0] = 1
    uncovered = np.arange(1, right, dtype=np.int32)
    index = np.arange(-1, right - 1, dtype=np.int32)
    scores = np.zeros(primes[-1], dtype=np.int16)
    touched = np.zeros(primes[-1], dtype=np.int32)
    count = right - 1
    for i, prime in enumerate(primes):
        count = change(covers, uncovered, index, count, prime, residues[i], 1)
    best = residues.copy()
    best_count = count
    best_log = 1e20
    log_primes = np.log(primes.astype(np.float64))
    order = np.arange(len(primes))
    for iteration in range(rounds):
        if iteration > 0:
            # Leave most of the cover in place while allowing several
            # interacting residues to change in the same move.
            np.random.shuffle(order)
            cut = max(1, len(order) // np.random.randint(12, 38))
            for i in order[:cut]:
                count = change(covers, uncovered, index, count, primes[i], residues[i], -1)
                residues[i] = -1
            for i in order[:cut]:
                prime = primes[i]
                r = best_residue(uncovered, count, prime, right, scores, touched)
                residues[i] = r
                count = change(covers, uncovered, index, count, prime, r, 1)
        for sweep in range(sweeps):
            if sweep or iteration:
                np.random.shuffle(order)
            for i in order:
                prime = primes[i]
                count = change(covers, uncovered, index, count, prime, residues[i], -1)
                r = best_residue(uncovered, count, prime, right, scores, touched)
                residues[i] = r
                count = change(covers, uncovered, index, count, prime, r, 1)
        log_mod = log_primes[residues >= 0].sum()
        if count < best_count or (count == best_count and log_mod < best_log):
            best = residues.copy()
            best_count = count
            best_log = log_mod
        elif count > best_count + max(4, best_count // 40):
            for i, prime in enumerate(primes):
                if residues[i] != best[i]:
                    count = change(covers, uncovered, index, count, prime, residues[i], -1)
                    count = change(covers, uncovered, index, count, prime, best[i], 1)
                    residues[i] = best[i]
    return best_count, best_log + np.log(2), best


def load_start(path, primes, half_gap):
    result = np.full(len(primes), -1, dtype=np.int64)
    if path.exists():
        old = json.loads(path.read_text())
        if old["gap"] == 2 * half_gap:
            classes = old["offsets"]
            for i, p in enumerate(primes):
                result[i] = classes.get(str(p), -1)
    return result


def save(path, primes, residues, half_gap):
    offsets = {int(p): int(r) for p, r in zip(primes, residues) if r >= 0}
    sites = np.ones(half_gap, dtype=np.bool_)
    sites[0] = False
    for p, r in offsets.items():
        if r in (0, half_gap % p):
            raise ArithmeticError("An endpoint was covered")
        sites[r::p] = False
    uncovered = (2 * np.flatnonzero(sites)).tolist()
    path.write_text(json.dumps({"gap": 2 * half_gap, "offsets": offsets, "uncovered": uncovered}))
    return len(uncovered)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gap", type=int, default=122290)
    parser.add_argument("--limits", type=int, nargs="+", default=[10000, 13000, 16000, 18000])
    parser.add_argument("--rounds", type=int, default=20)
    parser.add_argument("--batches", type=int, default=3)
    args = parser.parse_args()
    half_gap = args.gap // 2
    if args.gap < 4 or args.gap % 2:
        parser.error("gap must be a positive even integer of at least 4")
    previous = Path("work/nonexistent_mask.json")
    for limit in args.limits:
        path = Path(f"work/mask_{args.gap}_{limit}.json")
        primes = np.array(list(primerange(3, limit + 1)), dtype=np.int64)
        start = load_start(path if path.exists() else previous, primes, half_gap)
        winner = (half_gap, math.inf)
        if path.exists():
            saved = json.loads(path.read_text())
            winner = (len(saved["uncovered"]), sum(math.log(int(q)) for q in saved["offsets"]) + math.log(2))
        for batch in range(args.batches):
            begin = perf_counter()
            count, log_mod, found = optimize(primes, half_gap, start, args.rounds, 3, 853 + batch)
            print(f"{limit=} {batch=} {count=} ln W={log_mod:.1f} seconds={perf_counter()-begin:.2f}", flush=True)
            if (count, log_mod) < winner:
                if save(path, primes, found, half_gap) != count:
                    raise ArithmeticError("Recomputed cover does not match")
                winner = (count, log_mod)
            start = found
        previous = path


if __name__ == "__main__":
    main()
