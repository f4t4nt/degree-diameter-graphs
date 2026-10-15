"""Counting relaxation for general free-block architectures.

A free block is a subdivided 4-regular base graph in which every vertex is adjacent to or on every
edge, so every internal pair is within distance 4 using block edges alone: subdivided K5 (5 hubs,
10 ports, 15 vertices, 20 slots), subdivided octahedron K_{2,2,2} (6, 12, 18, 24) and subdivided
K4,4 (8, 16, 24, 32). The record architecture is three K4,4 blocks. This scans multisets of free
blocks and bounds, by the same hub-to-port identity and hub-to-hub coverage used for the K4,4 family,
the largest connector count the counting allows:

  hub a of block i must reach the 2 b_j ports of block j:  sum_c x_i(c) f_j(c) >= b_i b_j
  (the b_i hubs of block i each need 2 b_j reach units and each port of c in block i serves 2 hubs)
  hub pairs (i, j) must be covered by single connectors:    sum_c x_i(c) x_j(c) >= b_i b_j / 4
  slots: sum_c x_i(c) = 4 b_i;  degrees: connector edges = (4 - |x|) per connector

Connector profiles are integer vectors over the blocks with 1 <= |x| <= 4 (routers excluded here:
this is a discovery scan, not an impossibility proof). Output: for each block multiset, the
relaxation's maximum n = sum 3 b_i + number of connectors, with the optimal profile mix.
"""
import itertools
import json
import sys

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

BLOCK_B = {"K5": 5, "Oct": 6, "K44": 8}


def profiles(k):
    return [t for t in itertools.product(range(5), repeat=k) if 1 <= sum(t) <= 4]


def solve(blocks, time_limit=60, fixed_n=None):
    """Maximise the connector count (fixed_n None) or, for a fixed n, maximise the total budget slack
    sum over ordered block pairs of (sum_c x_i f_j - b_i b_j)."""
    b = [BLOCK_B[t] for t in blocks]
    k = len(b)
    TY = profiles(k)
    X = {t: sum(t) for t in TY}
    FREE = [t for t in TY if X[t] < 4]
    ET = [(t, u) for i, t in enumerate(FREE) for u in FREE[i:]]
    nv = len(TY) + len(ET)
    ci = {t: i for i, t in enumerate(TY)}
    ei = {e: len(TY) + i for i, e in enumerate(ET)}
    rows, lo, hi = [], [], []

    def add(r, l, h):
        rows.append(r), lo.append(l), hi.append(h)

    for j in range(k):
        r = np.zeros(nv)
        for t in TY:
            r[ci[t]] = t[j]
        add(r, 4 * b[j], 4 * b[j])
    for t in FREE:
        r = np.zeros(nv)
        r[ci[t]] = -(4 - X[t])
        for (a, c), idx in ei.items():
            if a == t and c == t:
                r[idx] += 2
            elif a == t or c == t:
                r[idx] += 1
        add(r, 0, 0)
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            r = np.zeros(nv)
            for t in TY:
                r[ci[t]] = t[i] * t[j]
            for (a, c), idx in ei.items():
                r[idx] = a[i] * c[j] + c[i] * a[j]
            add(r, b[i] * b[j], np.inf)
    for i in range(k):
        for j in range(i + 1, k):
            r = np.zeros(nv)
            for t in TY:
                r[ci[t]] = t[i] * t[j]
            add(r, (b[i] * b[j] + 3) // 4, np.inf)
    count = np.zeros(nv)
    for t in TY:
        count[ci[t]] = 1
    budget = np.zeros(nv)
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            for t in TY:
                budget[ci[t]] += t[i] * t[j]
            for (a, c), idx in ei.items():
                budget[idx] += a[i] * c[j] + c[i] * a[j]
    need = sum(b[i] * b[j] for i in range(k) for j in range(k) if i != j)
    if fixed_n is not None:
        add(count, fixed_n - 3 * sum(b), fixed_n - 3 * sum(b))
    obj = -(budget if fixed_n is not None else count)
    res = milp(obj, constraints=LinearConstraint(np.array(rows), lo, hi), integrality=np.ones(nv),
               bounds=Bounds(0, np.inf), options={"time_limit": time_limit})
    if not res.success:
        return {"blocks": blocks, "feasible": False, "status": int(res.status), "n": fixed_n}
    xs = np.rint(res.x).astype(int)
    return {"blocks": blocks, "feasible": True, "block_vertices": 3 * sum(b), "connectors": int(count @ xs),
            "max_n": 3 * sum(b) + int(count @ xs), "n": fixed_n,
            "budget": int(budget @ xs), "budget_needed": need, "slack": int(budget @ xs) - need,
            "counts": {str(t): int(xs[ci[t]]) for t in TY if xs[ci[t]]},
            "edges": {f"{a}-{c}": int(xs[idx]) for (a, c), idx in ei.items() if xs[idx]}}


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    out = []
    for k in range(2, kmax + 1):
        for combo in itertools.combinations_with_replacement(["K5", "Oct", "K44"], k):
            if 3 * sum(BLOCK_B[t] for t in combo) > 150:
                continue
            r = solve(list(combo))
            out.append(r)
            print(json.dumps({kk: r[kk] for kk in r if kk not in ("counts", "edges")}), file=sys.stderr, flush=True)
    json.dump(out, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
