"""Counting relaxation for the three-K4,4-block family (hub-to-port budget plus hub-to-hub coverage).

Connector types are port profiles (x1, x2, x3), 0 <= x1+x2+x3 <= 4, where x_j is the number of ports
of block j adjacent to the connector; the connector has 4 - sum(x) edges to other connectors.

For a diameter-4 member, every hub a of block i must reach all 16 ports of each other block j within
four steps, and the only such paths are a-p-c-q and a-p-c-c'-q, so
    S_j(a) = sum over a's 8 slot connectors c of [x_j(c) + sum_{c' ~ c} x_j(c')] >= 16.
Summing over the 16 hubs outside block j and over j gives the identity
    sum_{a,j} S_j(a) = sum_c self(c) + sum_{cc'} w(c,c') >= 768,
    self(c) = 2 (x^2 - sum_j x_j^2),  w(c,c') = 2 sum_j [(x - x_j) x'_j + (x' - x'_j) x_j].
Hub-to-hub coverage between blocks i and j needs at least 16 connectors carrying a port of each
(a-p-c-q-b is the only path shape), weighted by x_i * x_j.

The ILP maximises the connector count (n - 72) subject to these counting constraints, with the
edge-type counts as variables; it is a relaxation (no simplicity or arrangement constraints), so its
optimum is an upper bound on n for the family, and infeasibility at a given n is a proof.
"""
import itertools
import json
import sys

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

TYPES = [t for t in itertools.product(range(5), repeat=3) if sum(t) <= 4]
X = {t: sum(t) for t in TYPES}


def self_w(t):
    return 2 * (X[t] ** 2 - sum(v * v for v in t))


def edge_w(t, u):
    return 2 * sum((X[t] - t[j]) * u[j] + (X[u] - u[j]) * t[j] for j in range(3))


FREE = [t for t in TYPES if X[t] < 4]
EDGE_TYPES = [(t, u) for i, t in enumerate(FREE) for u in FREE[i:]]
MIDDLE = [t for t in FREE if 4 - X[t] >= 2]
# ordered 2-paths (t0, t1, t2): c0 ~ c1 ~ c2 with c0 != c2; stored once per unordered end pair
PATH_TYPES = [(t0, t1, t2) for t1 in MIDDLE for i, t0 in enumerate(FREE) for t2 in FREE[i:]]


def misses(t):
    return sum(1 for v in t if v == 0)


def solve(fixed_n=None, maximize_budget=False, hub_connector=True, class_counts=None,
          per_block_pair=True):
    nv = len(TYPES) + len(EDGE_TYPES) + (len(PATH_TYPES) if hub_connector else 0)
    ci = {t: i for i, t in enumerate(TYPES)}
    ei = {e: len(TYPES) + i for i, e in enumerate(EDGE_TYPES)}
    pi = {p: len(TYPES) + len(EDGE_TYPES) + i for i, p in enumerate(PATH_TYPES)} if hub_connector else {}
    rows, lo, hi = [], [], []

    def add(row, l, h):
        rows.append(row)
        lo.append(l)
        hi.append(h)

    def ordered_edge_coeff(row, a, b, coeff):
        """Add coeff times the number of ordered adjacent pairs (type a, type b)."""
        key = (a, b) if (a, b) in ei else (b, a)
        row[ei[key]] += coeff * (2 if a == b else 1)

    for j in range(3):
        r = np.zeros(nv)
        for t in TYPES:
            r[ci[t]] = t[j]
        add(r, 32, 32)
    for t in FREE:
        r = np.zeros(nv)
        r[ci[t]] = -(4 - X[t])
        for (a, b), k in ei.items():
            if a == t and b == t:
                r[k] += 2
            elif a == t or b == t:
                r[k] += 1
        add(r, 0, 0)
    budget = np.zeros(nv)
    for t in TYPES:
        budget[ci[t]] = self_w(t)
    for (a, b), k in ei.items():
        budget[k] = edge_w(a, b)
    add(budget, 768, np.inf)
    if per_block_pair:
        # the same identity restricted to the hubs of block i reaching the ports of block j:
        # sum_c x_i(c) f_j(c) >= 64 for every ordered pair (i, j)
        for i in range(3):
            for j in range(3):
                if i == j:
                    continue
                r = np.zeros(nv)
                for t in TYPES:
                    r[ci[t]] = t[i] * t[j]
                for (a, b), k in ei.items():
                    r[k] = a[i] * b[j] + b[i] * a[j]
                add(r, 64, np.inf)
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        r = np.zeros(nv)
        for t in TYPES:
            r[ci[t]] = t[i] * t[j]
        add(r, 16, np.inf)
    if class_counts is not None:
        for xval, cnt in class_counts.items():
            r = np.zeros(nv)
            for t in TYPES:
                if X[t] == xval:
                    r[ci[t]] = 1
            add(r, cnt, cnt)
    if hub_connector:
        # 2-path consistency: for each ordered edge type (t0, t1) with t1 a middle type,
        # the paths t0-t1-* number (d(t1)-1) times the ordered (t0, t1) adjacencies.
        for t1 in MIDDLE:
            for t0 in FREE:
                r = np.zeros(nv)
                for t2 in FREE:
                    key = (t0, t1, t2) if (t0, t1, t2) in pi else (t2, t1, t0)
                    r[pi[key]] += 1 if t0 != t2 else 2
                ordered_edge_coeff(r, t0, t1, -(4 - X[t1] - 1))
                add(r, 0, 0)
        # hub-to-connector demand for routers only: a router has no ports, so every one of the 24
        # hubs must reach it through a slot connector and at most two connector edges
        # (a-p-c0-r or a-p-c0-c1-r). Connectors with ports can also be reached through a shared
        # port (a-p-c0-q-c), which the type counts do not see, so they are not constrained here.
        router = (0, 0, 0)
        r = np.zeros(nv)
        r[ci[router]] -= 24
        for (a, b), k in ei.items():
            if a == router and b != router:
                r[k] += 2 * X[b]
            elif b == router and a != router:
                r[k] += 2 * X[a]
        for (t0, t1, t2), k in pi.items():
            if t2 == router and t0 != router:
                r[k] += 2 * X[t0]
            elif t0 == router and t2 != router:
                r[k] += 2 * X[t2]
        add(r, 0, np.inf)
    count = np.zeros(nv)
    for t in TYPES:
        count[ci[t]] = 1
    if fixed_n is not None:
        add(count, fixed_n - 72, fixed_n - 72)
    c = -(budget if maximize_budget else count)
    res = milp(c, constraints=LinearConstraint(np.array(rows), lo, hi),
               integrality=np.ones(nv), bounds=Bounds(0, np.inf),
               options={"time_limit": 600})
    if not res.success:
        return {"status": int(res.status), "message": res.message}
    xs = np.rint(res.x).astype(int)
    counts = {str(t): int(xs[ci[t]]) for t in TYPES if xs[ci[t]]}
    edges = {f"{a}-{b}": int(xs[k]) for (a, b), k in ei.items() if xs[k]}
    return {"n": 72 + int(count @ xs), "budget": int(budget @ xs), "counts": counts, "edges": edges}


def main():
    out = {"types": len(TYPES), "edge_types": len(EDGE_TYPES)}
    best = solve()
    out["max_n_relaxation"] = best
    per_n = {}
    for n in range(104, (best["n"] if best else 104) + 1):
        r = solve(fixed_n=n, maximize_budget=True)
        per_n[n] = None if r is None else {"max_budget": r["budget"], "slack": r["budget"] - 768,
                                          "counts": r["counts"], "edges": r["edges"]}
    out["per_n_max_budget"] = per_n
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
