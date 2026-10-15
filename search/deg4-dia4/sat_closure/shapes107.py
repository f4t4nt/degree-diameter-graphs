"""Shape enumeration for three K4,4 blocks with slack (n = 72 + N), from components107.py output.

Same constraints as shapes106.py (deficit sum 3N - 96 with quads counting -1,
potential sum >= 96 - 3N, gadget slots (g, g, g) with g = 32 - 2u, per ordered block pair rows
sum_c x_i f_j >= 64 with 4 per unit, port lemma in inequality form as a slot matching, and the
equality form for budget-tight shapes), with two changes:
  - quads: a shape may contain q quads when the potential allows (each costs >= 5); the deficit
    allowance in the signature recursion is total + (quads that still fit the potential);
  - pruning of the signature recursion by slot imbalance: the remaining gadgets have deficit R and
    cost budget P (quarter units); a gadget with deficit d and k triples has 2d + 3k slots and its
    triples cost at least 2 quarter units each, a triple pair has 6 slots for >= 8, so the remaining
    gadgets contribute at most 2R + 1.5P slots, which must cover 3 max(slots) - sum(slots).

Usage: shapes107.py N COMPONENTS.json > shapes.json
"""
import collections
import itertools
import json
import sys
import time
from fractions import Fraction as FR

import networkx as nx

TYPES_UNIT = {((0, 1, 2), (2, 1, 0)), ((0, 2, 1), (2, 0, 1)), ((1, 0, 2), (1, 2, 0)), ((1, 1, 1), (1, 1, 1))}
PAIRS = [(i, j) for i in range(3) for j in range(3) if i != j]


def is_unit(c):
    return c["n"] == 2 and tuple(sorted(tuple(t) for t in c["types"])) in TYPES_UNIT


def analyse(c):
    types = [tuple(t) for t in c["types"]]
    nbr = {i: set() for i in range(c["n"])}
    for a, b in c["edges"]:
        nbr[a].add(b)
        nbr[b].add(a)
    f = [tuple(types[i][j] + sum(types[k][j] for k in nbr[i]) for j in range(3)) for i in range(c["n"])]
    c["f"] = f
    c["deficit"] = sum(3 - sum(t) for t in types)
    c["pot4"] = int(FR(c["potential"]) * 4)
    c["rows"] = tuple(sum(types[h][i] * f[h][j] for h in range(c["n"])) for (i, j) in PAIRS)
    slots = []
    for h in range(c["n"]):
        for j in range(3):
            for _ in range(types[h][j]):
                need = tuple(4 - f[h][i] if i != j else 0 for i in range(3))
                slots.append((h, j, need))
    c["slots_list"] = slots
    c["sig"] = (c["deficit"], c["pot4"], tuple(c["slots"]), c["rows"], c["n"])
    return c


def matching_ok(gadgets, exact=False):
    """Port lemma as a matching problem on gadget slots; see shapes106.py for the derivation.
    exact: budget-tight shape, every (port, foreign block) pair has f_i(h1) + f_i(h2) = 4 exactly,
    co-holders are never adjacent, a common connector neighbour has no foreign port, two connectors
    co-hold at most one port; exhaustive search."""
    allslots = []
    for gi, g in enumerate(gadgets):
        nbr = {i: set() for i in range(len(g["types"]))}
        for a, b in g["edges"]:
            nbr[a].add(b)
            nbr[b].add(a)
        for (h, j, need) in g["slots_list"]:
            allslots.append(((gi, h), j, need, g["f"][h], nbr[h], [tuple(t) for t in g["types"]]))
    if exact:
        weak = [k for k, s in enumerate(allslots) if any(v != 2 for i, v in enumerate(s[2]) if i != s[1])]
    else:
        weak = [k for k, s in enumerate(allslots) if any(v > 2 for v in s[2])]
    if not weak:
        return True
    weakset = set(weak)

    def compatible(k, l):
        (c1, j1, need1, f1, n1, types), (c2, j2, need2, f2, n2, _) = allslots[k], allslots[l]
        if c1 == c2 or j1 != j2:
            return False
        if exact:
            if not (all(f2[i] == need1[i] for i in range(3) if i != j1) and all(f1[i] == need2[i] for i in range(3) if i != j2)):
                return False
            if c1[0] == c2[0]:
                if c2[1] in n1:
                    return False
                for c in n1 & n2:
                    if any(types[c][i] for i in range(3) if i != j1):
                        return False
            return True
        return all(f2[i] >= need1[i] for i in range(3) if i != j1) and all(f1[i] >= need2[i] for i in range(3) if i != j2)

    if exact:
        cands = {k: [l for l in range(len(allslots)) if l != k and compatible(k, l)] for k in weak}
        order = sorted(weak, key=lambda k: len(cands[k]))
        used, pairs_used = set(), set()

        def rec(idx):
            if idx == len(order):
                return True
            k = order[idx]
            if k in used:
                return rec(idx + 1)
            for l in cands[k]:
                if l in used:
                    continue
                key = frozenset((allslots[k][0], allslots[l][0]))
                if key in pairs_used:
                    continue
                used.update((k, l))
                pairs_used.add(key)
                if rec(idx + 1):
                    return True
                used.difference_update((k, l))
                pairs_used.discard(key)
            return False

        return rec(0)

    G = nx.Graph()
    for k in weak:
        G.add_node(k)
        for l in range(len(allslots)):
            if l != k and (l not in weakset or l > k) and compatible(k, l):
                G.add_edge(k, l, weight=2 if l in weakset else 1)
    if any(G.degree(k) == 0 for k in weak):
        return False
    M = nx.max_weight_matching(G, maxcardinality=False)
    covered = sum(1 for a, b in M for x in (a, b) if x in weakset)
    return covered == len(weak)


def main():
    t0 = time.time()
    N = int(sys.argv[1])
    data = json.load(open(sys.argv[2]))
    comps = data["components"] if isinstance(data, dict) else data
    bound4 = (96 - 3 * N) * 4
    total_def = 3 * N - 96
    gadgets = [analyse(c) for c in comps if not is_unit(c)]
    gadgets = [g for g in gadgets if g["pot4"] >= bound4]
    groups = collections.defaultdict(list)
    for g in gadgets:
        groups[g["sig"]].append(g)
    sigs = sorted(groups, key=lambda s: (s[0], -s[1]))
    print(f"{len(gadgets)} gadgets with potential >= {bound4}/4, {len(sigs)} signatures", file=sys.stderr, flush=True)
    last_quad = max([k for k, s in enumerate(sigs) if s[0] < 0], default=-1)
    stats = collections.Counter()
    sig_shapes = []

    def record(chosen, pot4, slots, conn, rows):
        rest = N - conn
        if rest % 2:
            stats["odd_rest"] += 1
            return
        u = rest // 2
        g = 32 - 2 * u
        if any(s != g for s in slots):
            stats["unbalanced"] += 1
            return
        stats["balanced"] += 1
        if any(4 * u + r < 64 for r in rows):
            stats["rows_fail"] += 1
            return
        stats["rows_ok"] += 1
        sig_shapes.append((u, pot4, list(chosen)))

    def rec(start, chosen, pot4, deff, slots, conn, rows):
        if deff == total_def and chosen:
            record(chosen, pot4, slots, conn, rows)
        P = pot4 - bound4  # remaining cost budget in quarter units
        allow = total_def + (P // 20 if start <= last_quad else 0)
        R = allow - deff
        need_fill = 3 * max(slots) - sum(slots)
        if need_fill > 2 * R + 1.5 * P:
            stats["pruned_imbalance"] += 1
            return
        for k in range(start, len(sigs)):
            s = sigs[k]
            if pot4 + s[1] < bound4:
                continue
            nd = deff + s[0]
            if nd > allow:
                continue
            if conn + s[4] > N:
                continue
            ns = [slots[0] + s[2][0], slots[1] + s[2][1], slots[2] + s[2][2]]
            if max(ns) > 32:
                continue
            rec(k, chosen + [s], pot4 + s[1], nd, ns, conn + s[4], [a + b for a, b in zip(rows, s[3])])

    rec(0, [], 0, 0, [0, 0, 0], 0, [0] * 6)
    print(json.dumps(dict(stats)), f"{time.time() - t0:.0f} s", file=sys.stderr, flush=True)
    shapes = []
    for u, pot4, chosen in sig_shapes:
        cnt = collections.Counter(chosen)
        pools = [list(itertools.combinations_with_replacement(groups[s], m)) for s, m in cnt.items()]
        for combo in itertools.product(*pools):
            gs = [g for part in combo for g in part]
            stats["gadget_multisets"] += 1
            if not matching_ok(gs):
                stats["matching_fail"] += 1
                continue
            if pot4 == bound4 and not matching_ok(gs, exact=True):
                stats["tight_exact_matching_fail"] += 1
                continue
            stats["survive"] += 1
            shapes.append({"units": u, "potential": str(FR(pot4, 4)), "budget": 768 + (pot4 - bound4) // 2,
                           "gadgets": [{"types": g["types"], "edges": g["edges"], "potential": g["potential"], "f": g["f"]} for g in gs]})
    print(json.dumps(dict(stats)), f"{time.time() - t0:.0f} s", file=sys.stderr, flush=True)
    shapes.sort(key=lambda s: (-FR(s["potential"]), len(s["gadgets"]), sum(len(g["types"]) for g in s["gadgets"])))
    json.dump({"N": N, "n": 72 + N, "bound": str(FR(bound4, 4)), "stats": dict(stats), "shapes": shapes}, sys.stdout, indent=0)


if __name__ == "__main__":
    main()
