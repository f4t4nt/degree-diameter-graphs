"""Which connector types and connector-edge types can occur in a 108-vertex Oct x4 graph.

Architecture: four subdivided octahedra (6 hubs, 12 ports each; hubs saturated; ports have two
free slots), 36 connectors, no port-port edges, maximum degree 4 (regularity is NOT assumed: a
connector type is (x, d) with x = ports per block, d = connector-connector degree, |x| + d <= 4,
and block slot sums are <= 24).

Hub-to-port identity (TEMPLATE_CALCULUS.md section 4, adapted): hub a of block i reaches a port q
of block j != i only along a-p-c-q or a-p-c-c'-q, so sum over the slot connectors c of a of
f_j(c) >= 12 where f_j(c) = x_j(c) + sum_{c' ~ c} x_j(c'). Summing over the 6 hubs of block i:
    row(i, j) = sum_c x_i(c) f_j(c) >= 36        (12 ordered block pairs)
and the sum of the 12 rows is  sum_c T(c), T(c) = |x||F| - x.F  (F = x + sum of neighbours' x),
which must be >= 432. Hub-to-hub coverage (a-p-c-q-b only): sum_c x_i x_j >= 9 per block pair.

This script maximises sum_c T(c) over type counts n_t and edge-type counts e_tu (integer),
then, for every variable, maximises again with that variable forced >= 1. A variable whose
forced maximum is below 432 cannot occur in any Oct x4 graph on 108 vertices. Output: the
surviving types and edge types (expected: straight and bent triples with d = 1, cross doubles
with d = 2, triple-triple edges with combined profile a permutation of (2,2,1,1), and
double-double edges between complementary cross types).
"""
import itertools
import json
import sys
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

NB = 4


def types():
    out = []
    for x in itertools.product(range(5), repeat=NB):
        if sum(x) > 4:
            continue
        for d in range(0, 5 - sum(x)):
            out.append((x, d))
    return out


def selfv(x):
    return sum(x) ** 2 - sum(v * v for v in x)


def v(x, u):
    return sum(x) * sum(u) - sum(a * b for a, b in zip(x, u))


def build():
    TY = types()
    FREE = [t for t in TY if t[1] >= 1]
    ET = [(FREE[i], FREE[j]) for i in range(len(FREE)) for j in range(i, len(FREE))]
    nv = len(TY) + len(ET)
    ti = {t: i for i, t in enumerate(TY)}
    ei = {e: len(TY) + i for i, e in enumerate(ET)}
    rows, lo, hi = [], [], []

    def add(r, l, h):
        rows.append(r)
        lo.append(l)
        hi.append(h)

    r = np.zeros(nv)
    for t in TY:
        r[ti[t]] = 1
    add(r, 36, 36)
    for j in range(NB):
        r = np.zeros(nv)
        for t in TY:
            r[ti[t]] = t[0][j]
        add(r, 0, 24)
    for t in FREE:
        r = np.zeros(nv)
        r[ti[t]] = -t[1]
        for (a, b), idx in ei.items():
            if a == t and b == t:
                r[idx] += 2
            elif a == t or b == t:
                r[idx] += 1
        add(r, 0, 0)
    obj = np.zeros(nv)
    for i in range(NB):
        for j in range(NB):
            if i == j:
                continue
            r = np.zeros(nv)
            for t in TY:
                r[ti[t]] = t[0][i] * t[0][j]
            for (a, b), idx in ei.items():
                r[idx] = a[0][i] * b[0][j] + b[0][i] * a[0][j]
            add(r, 36, np.inf)
            obj += r
    for i in range(NB):
        for j in range(i + 1, NB):
            r = np.zeros(nv)
            for t in TY:
                r[ti[t]] = t[0][i] * t[0][j]
            add(r, 9, np.inf)
    # sanity: the objective equals sum_t n_t self(t) + sum_e e 2 v(t,u)
    chk = np.zeros(nv)
    for t in TY:
        chk[ti[t]] = selfv(t[0])
    for (a, b), idx in ei.items():
        chk[idx] = 2 * v(a[0], b[0])
    assert np.array_equal(obj, chk)
    return TY, ET, ti, ei, np.array(rows), np.array(lo), np.array(hi), obj


def solve_milp(A, lo, hi, obj, lb):
    res = milp(-obj, constraints=LinearConstraint(A, lo, hi), integrality=np.ones(len(obj)),
               bounds=Bounds(lb, np.full(len(obj), np.inf)), options={"time_limit": 120})
    if res.status == 2:  # infeasible
        return None, None
    assert res.success, res
    return -res.fun, np.rint(res.x).astype(int)


def main():
    t0 = time.time()
    TY, ET, ti, ei, A, lo, hi, obj = build()
    nv = len(obj)
    names = [f"n{t}" for t in TY] + [f"e{a}-{b}" for a, b in ET]
    rec = {"types": len(TY), "edge_types": len(ET), "rows": int(A.shape[0])}
    lp = linprog(-obj, A_ub=np.vstack([A, -A]), b_ub=np.concatenate([hi, -lo]).clip(max=1e9),
                 bounds=(0, None), method="highs")
    assert lp.success
    rec["lp_max"] = round(-lp.fun, 6)
    red = lp.lower.marginals  # reduced costs for the minimisation; > 0 means forced to zero
    rec["lp_zero_reduced_cost"] = int(np.sum(red <= 1e-9))
    best, x = solve_milp(A, lo, hi, obj, np.zeros(nv))
    rec["ilp_max"] = best
    rec["ilp_argmax"] = {names[k]: int(x[k]) for k in range(nv) if x[k]}
    print(json.dumps(rec), file=sys.stderr, flush=True)
    assert abs(best - 432) < 1e-6
    # regularity: forcing a connector of degree < 4 or a missing slot
    # every variable is tested by an ILP (the LP reduced costs are reported only as a cross-check)
    lb = np.zeros(nv)
    survivors, excluded = [], []
    for k in range(nv):
        lbk = lb.copy()
        lbk[k] = 1
        val, _ = solve_milp(A, lo, hi, obj, lbk)
        if val is None or val < 432 - 1e-6:
            excluded.append((names[k], "infeasible" if val is None else round(float(val), 3), round(float(red[k]), 4)))
        else:
            survivors.append((names[k], round(float(red[k]), 4)))
    rec["survivors"] = survivors
    rec["n_excluded"] = len(excluded)
    rec["excluded"] = excluded
    rec["max_forced_value_among_excluded"] = max((e[1] for e in excluded if e[1] != "infeasible"), default=None)
    # slot deficiency: total slots <= 95
    r = np.zeros(nv)
    for t in TY:
        r[ti[t]] = sum(t[0])
    A2 = np.vstack([A, r])
    val, _ = solve_milp(A2, np.append(lo, 0), np.append(hi, 95), obj, np.zeros(nv))
    rec["max_with_slots_le_95"] = val
    rec["time_s"] = round(time.time() - t0, 1)
    print(json.dumps(rec), file=sys.stderr, flush=True)
    json.dump(rec, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
