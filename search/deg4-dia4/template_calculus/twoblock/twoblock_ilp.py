"""Counting relaxation for the two-block architecture (two vertex-disjoint subdivided K4,4 blocks,
saturated hubs, no port-port edge, arbitrary remainder W on m = n - 48 vertices).

W-vertex types are profiles t = (x1, x2), x_j = number of block-j ports adjacent to the vertex,
|t| <= 4, with 4 - |t| plain (W-W) edges. Variables: type counts n_t, plain-edge type counts e_{t,u},
2-path type counts p_{t0,t1,t2} (middle type with >= 2 plain edges). Constraints (all necessary,
RESULTS.md section 1):
  slots            sum_t t_j n_t = 32 (j = 1, 2)
  degrees          plain half-edges of type t = (4 - |t|) n_t
  2-path closure   paths t0-t1-* = (3 - |t1|) x ordered adjacencies (t0, t1)
  hub-to-port rows sum_c x_1 f_2 >= 64, sum_c x_2 f_1 >= 64            (B = 2 (row12 + row21) >= 256)
  hub-to-hub       sum_c x_1 x_2 >= 16
  router hubs      every router has >= 4 block-j port incidences within plain radius 2 (both j)
  pure singles     a (1,0) needs >= 1 block-2 port incidence within plain radius 2 beyond its
                   coholder's <= 3 (and symmetrically)
  hub-to-W         8 |Z_i| <= 2 sum_c x_i(c) g_i(c), g_i(c) <= sum_{c'~c} (4 - |x(c')|) + x_j(c)
  W-W Moore        sum_r |B''_4(r)| >= m^2 - 768 with the walk-count upper bound on the left
It is a relaxation (no arrangement constraints), so its optimum bounds n for the architecture, and
the max/min of B, cross incidences and routers at fixed n measure the slack of each constraint.

Usage: python twoblock_ilp.py  -> runs/twoblock_ilp.json
"""
import itertools
import json
import pathlib
import sys

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

HERE = pathlib.Path(__file__).resolve().parent
TYPES = [t for t in itertools.product(range(5), repeat=2) if sum(t) <= 4]
X = {t: sum(t) for t in TYPES}
FREE = [t for t in TYPES if X[t] < 4]
EDGE_TYPES = [(t, u) for i, t in enumerate(FREE) for u in FREE[i:]]
MIDDLE = [t for t in FREE if 4 - X[t] >= 2]
PATH_TYPES = [(t0, t1, t2) for t1 in MIDDLE for i, t0 in enumerate(FREE) for t2 in FREE[i:]]
ROUTER = (0, 0)


def solve(fixed_m=None, objective="m", sense="max", extra=None, time_limit=300):
    nv = len(TYPES) + len(EDGE_TYPES) + len(PATH_TYPES)
    ci = {t: i for i, t in enumerate(TYPES)}
    ei = {e: len(TYPES) + i for i, e in enumerate(EDGE_TYPES)}
    pi = {p: len(TYPES) + len(EDGE_TYPES) + i for i, p in enumerate(PATH_TYPES)}
    rows, lo, hi = [], [], []

    def add(row, l, h):
        rows.append(row)
        lo.append(l)
        hi.append(h)

    def ekey(a, b):
        return (a, b) if (a, b) in ei else (b, a)

    def ordered_adj(row, a, b, coeff):
        """coeff x (number of ordered adjacent pairs (type a, type b))."""
        row[ei[ekey(a, b)]] += coeff * (2 if a == b else 1)

    def pkey(t0, t1, t2):
        return (t0, t1, t2) if (t0, t1, t2) in pi else (t2, t1, t0)

    def ordered_path(row, t0, t1, t2, coeff):
        """coeff x (number of ordered 2-paths c0-c1-c2 with types t0, t1, t2, c0 != c2)."""
        row[pi[pkey(t0, t1, t2)]] += coeff * (2 if t0 == t2 else 1)

    # slots
    for j in range(2):
        r = np.zeros(nv)
        for t in TYPES:
            r[ci[t]] = t[j]
        add(r, 32, 32)
    # degrees
    for t in FREE:
        r = np.zeros(nv)
        r[ci[t]] = -(4 - X[t])
        for (a, b), k in ei.items():
            if a == t and b == t:
                r[k] += 2
            elif a == t or b == t:
                r[k] += 1
        add(r, 0, 0)
    # 2-path closure
    for t1 in MIDDLE:
        for t0 in FREE:
            r = np.zeros(nv)
            for t2 in FREE:
                ordered_path(r, t0, t1, t2, 1)
            ordered_adj(r, t0, t1, -(4 - X[t1] - 1))
            add(r, 0, 0)
    # hub-to-port rows and budget
    rowvec = {}
    for i, j in [(0, 1), (1, 0)]:
        r = np.zeros(nv)
        for t in TYPES:
            r[ci[t]] = t[i] * t[j]
        for (a, b), k in ei.items():
            r[k] = a[i] * b[j] + b[i] * a[j]
        rowvec[(i, j)] = r
        add(r, 64, np.inf)
    budget = 2 * (rowvec[(0, 1)] + rowvec[(1, 0)])
    # hub-to-hub
    cross = np.zeros(nv)
    for t in TYPES:
        cross[ci[t]] = t[0] * t[1]
    add(cross, 16, np.inf)
    # routers: >= 4 block-j port incidences within plain radius 2
    for j in range(2):
        r = np.zeros(nv)
        r[ci[ROUTER]] -= 4
        for u in FREE:
            if u != ROUTER:
                ordered_adj(r, ROUTER, u, u[j])
        for t1 in MIDDLE:
            for t2 in FREE:
                if t2[j]:
                    ordered_path(r, ROUTER, t1, t2, t2[j])
        add(r, 0, np.inf)
    # pure singles: block-(1-j) hubs via plain radius 2 or the coholder (<= 3 foreign ports)
    for j in range(2):
        t = (1, 0) if j == 0 else (0, 1)
        r = np.zeros(nv)
        r[ci[t]] -= 4 - 3
        for u in FREE:
            if u[1 - j]:
                ordered_adj(r, t, u, u[1 - j])
        for t1 in MIDDLE:
            for t2 in FREE:
                if t2[1 - j]:
                    ordered_path(r, t, t1, t2, t2[1 - j])
        add(r, 0, np.inf)
    # hub-to-W: 8 |Z_i| <= 2 sum_c x_i(c) g_i(c)
    for i in range(2):
        j = 1 - i
        r = np.zeros(nv)
        for t in TYPES:
            if t[i] == 0:
                r[ci[t]] -= 8  # member of Z_i
            r[ci[t]] += 2 * t[i] * t[j]  # coholders of block-j ports
        for (a, b), k in ei.items():
            # ordered adjacency (c of type a, c' of type b): a_i (4 - |b|), and the reverse
            r[k] += 2 * (a[i] * (4 - X[b]) + b[i] * (4 - X[a])) * (1 if a != b else 1)
            # for a == b the unordered count e_{a,a} represents 2 ordered pairs, each giving a_i (4-|a|):
            # 2 * (a_i (4-|a|) + a_i (4-|a|)) already equals 4 a_i (4-|a|) = 2 ordered x 2 a_i (4-|a|). ok.
        add(r, 0, np.inf)
    # W-W Moore: sum_r |B''_4(r)| <= m + W1 + W2 + W3 + W4 (non-backtracking walks, weights plain 1, port 2)
    # and >= m^2 - 768 for the architecture; linear only at fixed m.
    count = np.zeros(nv)
    for t in TYPES:
        count[ci[t]] = 1
    if fixed_m is not None:
        r = np.zeros(nv)
        for t in TYPES:
            d = 4 - X[t]
            # 1 (itself) + W1 plain + W2 plain-plain + port-edge (|t|) + W3 (plain,port)+(port,plain)
            # + W4 (plain,port,plain) bound 3 |t| d + port-port |t| (|t| - 1)
            r[ci[t]] += 1 + d + d * (d - 1) + X[t] + 2 * X[t] * d + 3 * X[t] * d + X[t] * (X[t] - 1)
        for (a, b), k in ei.items():
            da, db = 4 - X[a], 4 - X[b]
            # W3 plain^3: ordered walks through the edge as middle: 2 (da-1)(db-1)
            # W4 (plain,plain,port) + (port,plain,plain): 2 [(da-1)|b| + (db-1)|a|]
            r[k] += 2 * (da - 1) * (db - 1) + 2 * ((da - 1) * X[b] + (db - 1) * X[a])
        for (t0, t1, t2), k in pi.items():
            # W4 plain^4: ordered walks r-c0-c1-c2-r' with c0-c1-c2 the path: 2 (d0-1)(d2-1) per path
            d0, d2 = 4 - X[t0], 4 - X[t2]
            r[k] += 2 * (d0 - 1) * (d2 - 1)
        add(r, fixed_m * fixed_m - 768, np.inf)
        add(count, fixed_m, fixed_m)
    if extra:
        for row, l, h in extra(ci, ei, pi, nv, ordered_adj, ordered_path):
            add(row, l, h)
    objs = {"m": count, "B": budget, "cross": cross, "routers": np.eye(nv)[ci[ROUTER]],
            "singles": np.eye(nv)[ci[(1, 0)]] + np.eye(nv)[ci[(0, 1)]],
            "quads": sum(np.eye(nv)[ci[t]] for t in TYPES if X[t] == 4)}
    c = objs[objective] * (-1 if sense == "max" else 1)
    res = milp(c, constraints=LinearConstraint(np.array(rows), lo, hi), integrality=np.ones(nv),
               bounds=Bounds(0, np.inf), options={"time_limit": time_limit})
    if not res.success:
        return {"status": int(res.status), "message": res.message}
    xs = np.rint(res.x).astype(int)
    return {"m": int(count @ xs), "n": 48 + int(count @ xs), "B": int(budget @ xs), "cross": int(cross @ xs),
            "routers": int(xs[ci[ROUTER]]), "value": int(objs[objective] @ xs),
            "counts": {str(t): int(xs[ci[t]]) for t in TYPES if xs[ci[t]]},
            "edges": {f"{a}-{b}": int(xs[k]) for (a, b), k in ei.items() if xs[k]}}


def main():
    out = {"types": len(TYPES), "edge_types": len(EDGE_TYPES), "path_types": len(PATH_TYPES)}
    best = solve(objective="m", sense="max")
    out["max_m_relaxation_without_WW"] = best
    print("max m (no W-W row):", best.get("m"), "n =", best.get("n"), file=sys.stderr, flush=True)
    per = {}
    for m in [56, 57, 58, 59, 60, 61, 62, 64, 68, 72, 80]:
        rec = {}
        for obj, sense in [("B", "max"), ("B", "min"), ("cross", "min"), ("routers", "min"), ("routers", "max"),
                           ("singles", "max"), ("quads", "max")]:
            r = solve(fixed_m=m, objective=obj, sense=sense)
            rec[f"{sense}_{obj}"] = r
            print(m, sense, obj, "->", r.get("value", r.get("message")), file=sys.stderr, flush=True)
        per[48 + m] = rec
    out["per_n"] = per
    # largest m for which the relaxation with the W-W row is feasible (binary search up to 160)
    feas_m = None
    for m in range(56, 161):
        r = solve(fixed_m=m, objective="B", sense="max", time_limit=60)
        if "value" in r:
            feas_m = m
        else:
            break
    out["max_m_feasible_with_WW_row"] = feas_m
    print("max m with the W-W row:", feas_m, "n =", (48 + feas_m) if feas_m else None, file=sys.stderr, flush=True)
    json.dump(out, open(HERE / "runs/twoblock_ilp.json", "w"), indent=1)


if __name__ == "__main__":
    main()
