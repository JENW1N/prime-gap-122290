"""Repair an incomplete CRT cover with a sparse set of larger primes."""

import argparse
from collections import defaultdict
import heapq
import json
from math import log
from pathlib import Path
from random import Random
from time import perf_counter

from sympy import primerange


def repair(mask, lower_prime, upper_prime, seed):
    start = perf_counter()
    half = mask["gap"] // 2
    sites = [j // 2 for j in mask["uncovered"]]
    prime_list = list(primerange(lower_prime + 1, max(upper_prime, half) + 1))
    primes = [p for p in prime_list if p <= upper_prime and p < half]
    # For p > sqrt(half), a pair of offsets can share at most one
    # admissible repair prime. Build the groups from differences.
    if lower_prime ** 2 < half:
        raise ValueError("The lower cutoff must exceed sqrt(gap/2)")
    factor = [0] * half
    for q in primes:
        for multiple in range(q, half, q):
            factor[multiple] = q
    groups = defaultdict(set)
    for i, a in enumerate(sites):
        for j in range(i + 1, len(sites)):
            b = sites[j]
            q = factor[b - a]
            if q:
                residue = a % q
                if residue != 0 and residue != half % q:
                    members = groups[q, residue]
                    members.add(i)
                    members.add(j)
    keys = list(groups)
    vertices = list(groups.values())
    log_q = [log(q) for q, _ in keys]
    del groups
    neighbors = [[] for _ in sites]
    for edge, members in enumerate(vertices):
        for vertex in members:
            neighbors[vertex].append(edge)

    rng = Random(seed)
    sizes = list(map(len, vertices))
    noise = [rng.random() * 0.012 for _ in vertices]
    heap = [(-n / log_q[i] + noise[i], i, n) for i, n in enumerate(sizes)]
    heapq.heapify(heap)
    uncovered = set(range(len(sites)))
    assigned = {int(q): r for q, r in mask["offsets"].items()}
    initial_log = sum(log(q) for q in assigned) + log(2)
    while heap:
        _, edge, size = heapq.heappop(heap)
        prime, residue = keys[edge]
        if size != sizes[edge] or prime in assigned:
            continue
        assigned[prime] = residue
        for vertex in vertices[edge]:
            if vertex in uncovered:
                uncovered.remove(vertex)
                for incident in neighbors[vertex]:
                    sizes[incident] -= 1
                    n = sizes[incident]
                    if n >= 2 and keys[incident][0] not in assigned:
                        heapq.heappush(heap, (-n / log_q[incident] + noise[incident], incident, n))

    # A fresh prime greater than the offsets always covers one site.
    # Use the smallest available prime that preserves the endpoints.
    if len(prime_list) < 2 * len(sites) + 10:
        prime_list = list(primerange(lower_prime + 1, max(upper_prime, half) + 2 * len(sites) ** 2 + 1000))
    for vertex in sorted(uncovered):
        site = sites[vertex]
        for prime in prime_list:
            if prime not in assigned and site % prime not in (0, half % prime):
                assigned[prime] = site % prime
                break
        else:
            raise ArithmeticError("Ran out of fresh repair primes")
    final_log = sum(log(q) for q in assigned) + log(2)
    print(f"seed {seed}: {len(keys)} available groups, {len(uncovered)} singleton "
          f"repairs, ln W {initial_log:.2f} -> {final_log:.2f} "
          f"in {perf_counter()-start:.2f}s", flush=True)
    return final_log, {"gap": mask["gap"], "offsets": assigned, "uncovered": []}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mask", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--limit", type=int, default=61145)
    parser.add_argument("--tries", type=int, default=3)
    args = parser.parse_args()
    mask = json.loads(args.mask.read_text())
    lower = max(map(int, mask["offsets"]))
    best = float("inf")
    for seed in range(args.tries):
        logarithm, result = repair(mask, lower, args.limit, seed + 853)
        if logarithm < best:
            best = logarithm
            args.output.write_text(json.dumps(result))


if __name__ == "__main__":
    main()
