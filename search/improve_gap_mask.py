"""Improve a complete CRT cover by solving small weighted covering problems."""

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import random
from time import perf_counter

import numpy as np
from ortools.sat.python import cp_model
from sympy import primerange


def uncovered_sites(half_gap, classes):
    uncovered = np.ones(half_gap, dtype=np.bool_)
    uncovered[0] = False
    for prime, residue in classes.items():
        if not 0 < residue < prime or residue == half_gap % prime:
            raise ValueError(f"An endpoint is blocked modulo {prime}")
        uncovered[residue::prime] = False
    return np.flatnonzero(uncovered).tolist()


def group_classes(sites, primes, factors, half_gap):
    groups = defaultdict(set)
    available = set(primes)
    for i, a in enumerate(sites):
        for j in range(i + 1, len(sites)):
            prime = factors[sites[j] - a]
            if prime and prime in available:
                residue = a % prime
                if residue and residue != half_gap % prime:
                    members = groups[prime, residue]
                    members.add(i)
                    members.add(j)

    for prime in primes:
        if prime * prime >= half_gap:
            break
        residues = defaultdict(list)
        forbidden = half_gap % prime
        for index, site in enumerate(sites):
            residue = site % prime
            if residue and residue != forbidden:
                residues[residue].append(index)
        for residue, indices in residues.items():
            if len(indices) >= 2:
                groups[prime, residue] = indices
    return groups


def solve(sites, primes, factors, half_gap, dropped, timeout, workers):
    groups = group_classes(sites, primes, factors, half_gap)
    model = cp_model.CpModel()
    by_prime = defaultdict(list)
    by_site = [[] for _ in sites]
    cover_options = {}
    group_covered = set()

    for (prime, residue), indices in groups.items():
        take = model.new_bool_var("")
        by_prime[prime].append(take)
        cover_options[prime, residue] = take
        for index in indices:
            by_site[index].append(take)
        selected = dropped.get(prime) == residue
        if selected:
            group_covered.update(indices)
        model.add_hint(take, selected)

    remaining = [index for index in range(len(sites)) if index not in group_covered]
    singleton_hints = set()
    for prime, residue in dropped.items():
        if (prime, residue) not in cover_options:
            for index in remaining:
                if sites[index] % prime == residue:
                    singleton_hints.add(prime)
                    break

    reserve = {}
    objective = []
    for prime in primes:
        single = model.new_bool_var("")
        reserve[prime] = single
        model.add_hint(single, prime in singleton_hints)
        choices = by_prime[prime] + [single]
        model.add_at_most_one(choices)
        objective.append(round(10_000 * math.log(prime)) * sum(choices))

    # Each remaining site requires a singleton slot. Assign those slots by
    # matching afterwards, where each site forbids only its endpoint factors.
    slots = []
    for index, choices in enumerate(by_site):
        single = model.new_bool_var("")
        model.add_hint(single, index not in group_covered)
        model.add_bool_or(choices + [single])
        slots.append(single)
    model.add(sum(slots) <= sum(reserve.values()))
    model.minimize(sum(objective))

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = timeout
    solver.parameters.num_search_workers = workers
    status = solver.solve(model)
    if status not in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        return None, len(groups)

    classes = {
        (prime): residue for (prime, residue), take in cover_options.items() if solver.value(take)
    }
    left = [site for site in sites if all(site % q != r for q, r in classes.items())]
    spares = [prime for prime, var in reserve.items() if solver.value(var)]
    assigned = {}

    def match(site, tried):
        for prime in spares:
            if prime in tried or site % prime in (0, half_gap % prime):
                continue
            tried.add(prime)
            if prime not in assigned or match(assigned[prime], tried):
                assigned[prime] = site
                return True
        return False

    for site in left:
        if not match(site, set()):
            for prime in primes:
                if prime not in classes and prime not in spares:
                    spares.append(prime)
                    if match(site, set()):
                        break
            else:
                return None, len(groups)
    classes.update({prime: site % prime for prime, site in assigned.items()})
    return classes, len(groups)


def cost(classes):
    return math.log(2) + sum(map(math.log, classes))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--rounds", type=int, default=20)
    parser.add_argument("--timeout", type=float, default=6)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=853)
    args = parser.parse_args()

    original = json.loads(args.source.read_text())
    half_gap = original["gap"] // 2
    best = {int(q): int(r) for q, r in original["offsets"].items()}
    if original["uncovered"] or uncovered_sites(half_gap, best):
        raise ValueError("Start from a complete cover")

    rng = random.Random(args.seed)
    primes = list(primerange(3, max(36000, max(best) + 100)))
    factors = [0] * half_gap
    for prime in primes:
        if prime * prime > half_gap:
            for k in range(prime, half_gap, prime):
                factors[k] = prime

    begin = perf_counter()
    best_cost = cost(best)
    for turn in range(args.rounds):
        minimum = rng.choice([200, 450, 2000, 6000, 10000])
        removable = [q for q in best if q >= minimum]
        width = rng.randrange(75, 310)
        drop = set(rng.sample(removable, min(width, len(removable))))
        # Make the expensive tail available in every neighborhood.
        drop.update(sorted(best, reverse=True)[: rng.randrange(30, 70)])
        fixed = {q: r for q, r in best.items() if q not in drop}
        sites = uncovered_sites(half_gap, fixed)
        if not sites or len(sites) > 950:
            continue

        available = [q for q in primes if q not in fixed]
        chosen, count = solve(sites, available, factors, half_gap, {q: best[q] for q in drop}, args.timeout, args.workers)
        if chosen is None:
            print(f"{turn=} sites={len(sites)} options={count} no incumbent", flush=True)
            continue
        candidate = fixed | chosen
        value = cost(candidate)
        print(f"{turn=} sites={len(sites)} options={count} ln W={value:.2f} best={best_cost:.2f}", flush=True)
        if value < best_cost - 0.00001:
            if uncovered_sites(half_gap, candidate):
                raise ArithmeticError("Improved cover is incomplete")
            best, best_cost = candidate, value
            result = {"gap": 2 * half_gap, "offsets": best, "uncovered": []}
            args.output.write_text(json.dumps(result))
    print(f"Complete cover ln W={best_cost:.2f} after {perf_counter()-begin:.2f}s", flush=True)


if __name__ == "__main__":
    main()
