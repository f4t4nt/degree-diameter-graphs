"""Shape enumeration for three K4,4 blocks with slack: n = 72 + N, N connectors on 96 slots.

A shape is a multiset of gadgets (connected typed connector components that are not complementary
triple pairs, from potential3.py) plus u units (complementary triple pairs, two ports per block,
both members with carried vector f = (2,2,2)). Constraints:
  - deficit: sum over gadget connectors of (3 - |x|) = 3N - 96          (units have deficit 0)
  - potential: sum of gadget potentials >= 96 - 3N                      (units have potential 0)
  - slots: gadget slots are (g, g, g) with g = 32 - 2u, 2u + gadget connectors = N
  - per ordered block pair (i, j): sum_c x_i(c) f_j(c) >= 64, a unit contributes 4
  - port lemma (inequality form, TEMPLATE_CALCULUS.md section 7 / 12): at a port q of block j held
    by h1, h2 and for every block i != j, f_i(h1) + f_i(h2) >= 4 (2 (f_i(h1) + f_i(h2)) routes must
    cover the 8 hubs of block i). A unit member has f_i = 2, so a gadget slot (h, j) with f_i(h) <= 1
    for some i != j must be co-held by a gadget connector h2 with a block-j port and
    f_i(h2) >= 4 - f_i(h) for every such i; the weak slots must all be matched to distinct
    compatible slots of other connectors in the same block (a port has exactly two holders): a
    maximum-weight matching problem on the slot compatibility graph (networkx).
  - hub-to-hub coverage is left to the SAT model (it depends on the unit types).
Gadgets are grouped by coarse signature (deficit, potential, slots, rows, size); the enumeration
runs over signature multisets, and surviving signature multisets are expanded to gadget multisets
for the matching check. Output: the surviving shapes with their gadgets, for typed_slack_sat.py.

Usage: shapes106.py N COMPONENTS.json > shapes.json
"""
import collections
import itertools
import json
import sys
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
    c["nweak"] = sum(1 for s in slots if any(v > 2 for v in s[2]))
    c["sig"] = (c["deficit"], c["pot4"], tuple(c["slots"]), c["rows"], c["n"])
    return c


def matching_ok(gadgets, exact=False):
    """exact: budget-tight shape (total potential = bound), where every (port, foreign block) pair
    has f_i(h1) + f_i(h2) = 4 exactly (the total excess is an identity in the shape), so every slot
    with f_i(h) != 2 for some foreign i needs a co-holder with f_i(h2) = 4 - f_i(h) exactly, and the
    four carried block-i ports of N[h1] + N[h2] are distinct (a perfect matching): hence (a) h1, h2
    are not adjacent (each would carry its own foreign ports twice), (b) a common connector
    neighbour has no foreign port, (c) two connectors co-hold at most one port (a second shared
    port would be carried twice at the first). In exact mode the search is exhaustive with (c)."""
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
    N = int(sys.argv[1])
    comps = json.load(open(sys.argv[2]))
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
        allow = total_def + (1 if start <= last_quad else 0)
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
    print(json.dumps(dict(stats)), file=sys.stderr, flush=True)
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
    print(json.dumps(dict(stats)), file=sys.stderr, flush=True)
    shapes.sort(key=lambda s: (-FR(s["potential"]), len(s["gadgets"]), sum(len(g["types"]) for g in s["gadgets"])))
    json.dump({"N": N, "n": 72 + N, "bound": str(FR(bound4, 4)), "stats": dict(stats), "shapes": shapes}, sys.stdout, indent=0)


if __name__ == "__main__":
    main()
