"""Re-run the shape enumeration of shapes_fast.py with assumption (c) of the tight
equality matching removed.

matching_ok(exact=True) in shapes107.py (and shapes106.py) assumes
  (a) co-holders of a port are never adjacent,
  (b) a common connector neighbour of two co-holders has no foreign port,
  (c) two connectors co-hold at most one port.
(a) and (b) follow from hub-to-port exactness (a repeated member of N[h1] + N[h2] carries its
foreign ports twice). (c) does not: the second shared port lies in the same block j as the first,
and the port lemma at a block-j port only constrains foreign ports. So (c) may have excluded
genuine tight shapes at 106 and 107. This driver monkeypatches matching_ok to drop the
`pairs_used` test and reruns the (otherwise unchanged) pipeline; the shape lists are compared
with the originals by the caller (compare_shapes.py).

Usage: recheck_tight.py N --out=runs/shapes<n>_noc.json [--log=FILE] [--nosym]
"""
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import shapes107 as S  # noqa: E402
import shapes_fast  # noqa: E402

import networkx as nx  # noqa: E402


def matching_ok_noc(gadgets, exact=False):
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
        used = set()

        def rec(idx):
            if idx == len(order):
                return True
            k = order[idx]
            if k in used:
                return rec(idx + 1)
            for l in cands[k]:
                if l in used:
                    continue
                used.update((k, l))
                if rec(idx + 1):
                    return True
                used.difference_update((k, l))
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


S.matching_ok = matching_ok_noc

if __name__ == "__main__":
    shapes_fast.main()
