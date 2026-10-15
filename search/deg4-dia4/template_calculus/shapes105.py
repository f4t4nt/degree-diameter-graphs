"""Exhaustive counting analysis for n = 105 in the three-K4,4-block family (no port-port edges).

A 105-vertex member has 33 connectors on 96 port slots, so sum_c (3 - x(c)) = 3: with Q quads,
T triples, D doubles, S singles and R routers, D + 2S + 3R - Q = 3.

A "normal pair" is two triples joined by a connector edge whose port multiset has exactly two ports
in every block; both members then have reach vector f = (2,2,2). Every other connector is special:
quads, doubles, singles, routers, and the triples attached to them ("stubs"). For each shape, every
labelled configuration of specials is enumerated (profiles, connector edges among the non-triple
specials, stub-profile multisets per owner) and tested with sound necessary conditions:

1. hub-to-port budget per ordered block pair (i, j): sum_c x_i(c) f_j(c) >= 64, where
   f_j(c) = x_j(c) + sum_{c' ~ c} x_j(c'); a normal pair contributes exactly 4;
2. per-hub compensation: a hub holds 8 slot connectors and needs sum f_j >= 16 for both other
   blocks; normal-pair triples contribute exactly 2 each, so the special incidences of a block (with
   deficits f - 2) must be placeable in its 4x4 port grid, two per port, distinct ports per connector,
   with all row and column sums >= 0 in both coordinates;
3. hub-to-connector: a connector missing block i needs >= 4 ports of block i on connectors adjacent
   to it, at connector-distance two, or sharing one of its ports (unknown sharers bounded by 2);
4. hub-to-hub: >= 16 connectors carrying a port of both blocks, per block pair.

Shapes are first screened with the counting relaxation (counting_ilp.solve with class counts).
A configuration rejected here cannot occur in a diameter-4 member; nothing searches over graphs.
"""
import collections
import functools
import itertools
import json
import sys

import counting_ilp as C

TYPES = C.TYPES
X = C.X
TRIPLES = [t for t in TYPES if X[t] == 3]
QUADS = [t for t in TYPES if X[t] == 4]
DEG = {"D": 2, "S": 3, "R": 4}
PROFILES = {"D": [t for t in TYPES if X[t] == 2], "S": [t for t in TYPES if X[t] == 1], "R": [(0, 0, 0)]}
PAIRS = [(i, j) for i in range(3) for j in range(3) if i != j]


def shapes():
    out = []
    for Q in range(0, 34):
        for R in range(0, 2):
            for S in range(0, 4):
                D = 3 + Q - 2 * S - 3 * R
                T = 33 - Q - D - S - R
                if D >= 0 and T >= 0:
                    out.append((Q, T, D, S, R))
    return out


def skeletons(Q, T, D, S, R):
    """Non-triple specials with profiles, the simple graph among them, and stub counts."""
    units = ["D"] * D + ["S"] * S + ["R"] * R
    m = len(units)
    pairs = list(itertools.combinations(range(m), 2))
    for profiles in itertools.product(*[PROFILES[u] for u in units]):
        for mask in range(1 << len(pairs)):
            edges = [pairs[k] for k in range(len(pairs)) if mask >> k & 1]
            d = [0] * m
            for a, b in edges:
                d[a] += 1
                d[b] += 1
            if any(d[i] > DEG[units[i]] for i in range(m)):
                continue
            stubs = [DEG[units[i]] - d[i] for i in range(m)]
            if sum(stubs) > T or (T - sum(stubs)) % 2:
                continue
            yield units, profiles, edges, stubs


@functools.lru_cache(maxsize=None)
def placement_ok(incidences):
    """incidences: canonical tuple of (connector id, dj, dk). Can they be placed in the 4x4 port
    grid (two slots per port, distinct ports per connector) with all row and column deficit sums
    >= 0 in both coordinates? Solved as a small binary feasibility program (HiGHS)."""
    if not incidences:
        return True
    if sum(d[1] for d in incidences) < 0 or sum(d[2] for d in incidences) < 0:
        return False
    if all(d[1] >= 0 and d[2] >= 0 for d in incidences):
        return True
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp

    cells = [(r, c, s) for r in range(4) for c in range(4) for s in range(2)]
    m, nc = len(incidences), len(cells)
    nv = m * nc
    rows, lo, hi = [], [], []

    def add(row, l, h):
        rows.append(row), lo.append(l), hi.append(h)

    for a in range(m):
        r = np.zeros(nv)
        r[a * nc:(a + 1) * nc] = 1
        add(r, 1, 1)
    for k in range(nc):
        r = np.zeros(nv)
        for a in range(m):
            r[a * nc + k] = 1
        add(r, 0, 1)
    by_id = collections.defaultdict(list)
    for a, (cid, _, _) in enumerate(incidences):
        by_id[cid].append(a)
    for cid, members in by_id.items():
        if len(members) < 2:
            continue
        for pr in range(4):
            for pc in range(4):
                r = np.zeros(nv)
                for a in members:
                    for k, (cr, cc, _) in enumerate(cells):
                        if cr == pr and cc == pc:
                            r[a * nc + k] = 1
                add(r, 0, 1)
    for coord in (1, 2):
        for line in range(4):
            for which in (0, 1):  # rows then columns
                r = np.zeros(nv)
                for a, inc in enumerate(incidences):
                    for k, cell in enumerate(cells):
                        if cell[which] == line:
                            r[a * nc + k] = inc[coord]
                add(r, 0, np.inf)
    res = milp(np.zeros(nv), constraints=LinearConstraint(np.array(rows), lo, hi),
               integrality=np.ones(nv), bounds=Bounds(0, 1))
    return bool(res.success)


def canonical(inc):
    ids = {}
    out = []
    for cid, dj, dk in sorted(inc, key=lambda t: (t[1], t[2], str(t[0]))):
        ids.setdefault(cid, len(ids))
        out.append((ids[cid], dj, dk))
    return tuple(out)


def owner_options(o, special_nbr_sum, n_stubs):
    """For owner profile o (tuple) with fixed special-neighbour x-sum, enumerate stub multisets.
    Returns list of (contrib dict over PAIRS, f_owner, [f_stub...], stub profiles)."""
    opts = []
    for combo in itertools.combinations_with_replacement(TRIPLES, n_stubs):
        f_o = [o[j] + special_nbr_sum[j] + sum(t[j] for t in combo) for j in range(3)]
        contrib = {}
        for i, j in PAIRS:
            v = o[i] * f_o[j]
            for t in combo:
                v += t[i] * (t[j] + o[j])
            contrib[(i, j)] = v
        f_stubs = [[t[j] + o[j] for j in range(3)] for t in combo]
        opts.append((contrib, f_o, f_stubs, combo))
    return opts


def evaluate_full(units, profiles, adj, stub_choice, quads, normal_pairs):
    """Build the full special set and run checks 2-4 (check 1 already passed)."""
    prof = {i: tuple(profiles[i]) for i in range(len(units))}
    cadj = {i: set(adj[i]) for i in range(len(units))}
    nxt = len(units)
    for owner, (_, _, _, combo) in zip(range(len(units)), stub_choice):
        for t in combo:
            prof[nxt] = t
            cadj[nxt] = {owner}
            cadj[owner].add(nxt)
            nxt += 1
    for q in quads:
        prof[nxt] = q
        cadj[nxt] = set()
        nxt += 1
    f = {c: [prof[c][j] + sum(prof[c2][j] for c2 in cadj[c]) for j in range(3)] for c in prof}
    for i in range(3):
        j, k = [b for b in range(3) if b != i]
        inc = []
        for c in prof:
            for _ in range(prof[c][i]):
                inc.append((c, f[c][j] - 2, f[c][k] - 2))
        if not placement_ok(canonical(inc)):
            return "per-hub"
    for c in prof:
        for i in range(3):
            if prof[c][i]:
                continue
            ball = set(cadj[c])
            for c2 in list(ball):
                ball |= cadj[c2]
            ball.discard(c)
            if sum(prof[c2][i] for c2 in ball) + 2 * X[prof[c]] < 4:
                return "hub-connector"
    for i, j in ((0, 1), (0, 2), (1, 2)):
        if 2 * normal_pairs + sum(prof[c][i] * prof[c][j] for c in prof) < 16:
            return "hub-hub"
    return "pass"


def main():
    out = {"shapes": [], "survivors": [], "relaxation_infeasible_shapes": []}
    for Q, T, D, S, R in shapes():
        feas = C.solve(class_counts={4: Q, 3: T, 2: D, 1: S, 0: R}, fixed_n=105, maximize_budget=True)
        if "n" not in feas:
            out["relaxation_infeasible_shapes"].append([Q, T, D, S, R])
            continue
        rec = {"Q": Q, "T": T, "D": D, "S": S, "R": R, "relaxation_max_budget": feas["budget"],
               "skeletons": 0, "configurations_budget_checked": 0, "outcomes": collections.Counter()}
        for units, profiles, edges, stubs in skeletons(Q, T, D, S, R):
            rec["skeletons"] += 1
            m = len(units)
            adj = {i: set() for i in range(m)}
            for a, b in edges:
                adj[a].add(b)
                adj[b].add(a)
            normal_pairs = (T - sum(stubs)) // 2
            nbr_sum = [[sum(profiles[c2][j] for c2 in adj[i]) for j in range(3)] for i in range(m)]
            options = [owner_options(tuple(profiles[i]), nbr_sum[i], stubs[i]) for i in range(m)]
            ub = [{p: max(o[0][p] for o in options[i]) for p in PAIRS} for i in range(m)]
            suffix_ub = [{p: sum(ub[k][p] for k in range(i, m)) for p in PAIRS} for i in range(m + 1)]
            for quads in itertools.combinations_with_replacement(QUADS, Q):
                base = {p: 4 * normal_pairs + sum(q[p[0]] * q[p[1]] for q in quads) for p in PAIRS}

                def rec_owner(i, acc, chosen):
                    if any(acc[p] + suffix_ub[i][p] < 64 for p in PAIRS):
                        return
                    if i == m:
                        rec["configurations_budget_checked"] += 1
                        verdict = evaluate_full(units, profiles, adj, chosen, quads, normal_pairs)
                        rec["outcomes"][verdict] += 1
                        if verdict == "pass":
                            out["survivors"].append({
                                "shape": [Q, T, D, S, R], "units": units, "profiles": list(profiles),
                                "edges": edges, "stubs": [list(c[3]) for c in chosen], "quads": list(quads)})
                        return
                    for opt in options[i]:
                        acc2 = {p: acc[p] + opt[0][p] for p in PAIRS}
                        chosen.append(opt)
                        rec_owner(i + 1, acc2, chosen)
                        chosen.pop()

                rec_owner(0, dict(base), [])
        rec["outcomes"] = dict(rec["outcomes"])
        out["shapes"].append(rec)
        print(f"shape Q={Q} T={T} D={D} S={S} R={R}: relax {feas['budget']}, skeletons {rec['skeletons']}, "
              f"budget-feasible configs {rec['configurations_budget_checked']}, outcomes {rec['outcomes']}",
              file=sys.stderr, flush=True)
    out["survivor_count"] = len(out["survivors"])
    print(json.dumps(out, indent=1, default=list))


if __name__ == "__main__":
    main()
