"""Corrected enumeration of hub-to-port-tight unit designs (skeletons).

The first version (unit_designs.py / design_canon.py) treated every automorphism of the label
hierarchy (cycles, sides, labels) as a symmetry. Only the K4,4-induced label group is a symmetry; it
has order 64, not 128, and does not contain single within-side label swaps. That made the dedup too
coarse and the normalisation of labels lossy, which validate/classify104 exposed (26 known classes
were not splittings of the chosen representative). This script redoes the enumeration with the true
group and reduces designs by nauty certificates of the full unit-port-hub incidence graph.
"""
import collections
import itertools
import json
import resource
import sys

import pynauty

resource.setrlimit(resource.RLIMIT_AS, (6 * 1024 ** 3, 6 * 1024 ** 3))

N = 16
LABELS = ["A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2"]
POS = {l: i for i, l in enumerate(LABELS)}
SIDE = {"A1": "A", "A2": "A", "B1": "B", "B2": "B", "C1": "C", "C2": "C", "D1": "D", "D2": "D"}
EDGES = {
    "A1": ((0, 0), (1, 1)), "A2": ((0, 1), (1, 0)),
    "B1": ((2, 2), (3, 3)), "B2": ((2, 3), (3, 2)),
    "C1": ((0, 2), (1, 3)), "C2": ((0, 3), (1, 2)),
    "D1": ((2, 0), (3, 1)), "D2": ((2, 1), (3, 0)),
}
SIDES1 = {"A": [0, 1, 2, 3], "B": [4, 5, 6, 7], "C": [8, 9, 10, 11], "D": [12, 13, 14, 15]}
TAU1 = tuple(u ^ 1 for u in range(N))
OPP = {"A": "B", "B": "A", "C": "D", "D": "C"}


def label_group():
    inv = {frozenset(v): k for k, v in EDGES.items()}
    perms = set()
    for pl in itertools.permutations(range(4)):
        for pr in itertools.permutations(range(4)):
            for swap in (False, True):
                img = []
                for l in LABELS:
                    es = frozenset((pr[r], pl[a]) if swap else (pl[a], pr[r]) for a, r in EDGES[l])
                    if es not in inv:
                        break
                    img.append(inv[es])
                else:
                    perms.add(tuple(img))
    return sorted(perms)


def compose(p, q):
    return tuple(p[q[x]] for x in range(N))


def inverse(p):
    inv = [0] * N
    for i, v in enumerate(p):
        inv[v] = i
    return tuple(inv)


def conj(g, t):
    gi = inverse(g)
    return tuple(g[t[gi[x]]] for x in range(N))


def block1_stabiliser_true(L):
    """Unit permutations induced by the true label group on block 1's canonical layout, combined
    with swapping the two units of any twin pair."""
    out = set()
    for lam in L:
        for s in range(256):
            g = [0] * N
            for k in range(8):
                nk = POS[lam[k]]
                for eps in (0, 1):
                    g[2 * k + eps] = 2 * nk + (eps ^ ((s >> k) & 1))
            out.add(tuple(g))
    return sorted(out)


def valid_taus():
    out = []
    for pa in itertools.permutations(SIDES1["B"]):
        for pc in itertools.permutations(SIDES1["D"]):
            t = list(range(N))
            for a, b in zip(SIDES1["A"], pa):
                t[a], t[b] = b, a
            for c, d in zip(SIDES1["C"], pc):
                t[c], t[d] = d, c
            out.append(tuple(t))
    return out


def orbits_bipartite(t1, t2):
    seen = set()
    res = []
    for u in range(N):
        if u in seen:
            continue
        comp = {u: 0}
        stack = [u]
        while stack:
            x = stack.pop()
            for t in (t1, t2):
                y = t[x]
                if y not in comp:
                    comp[y] = 1 - comp[x]
                    stack.append(y)
        seen |= set(comp)
        res.append(comp)
    return res


def side_structures(tj, tk, ti):
    """All side maps for a block (no normalisation): orbits of <tj,tk> split 8+8 into Z1/Z2, flips
    per orbit, ti preserves sides."""
    orbs = orbits_bipartite(tj, tk)
    k = len(orbs)
    out = []
    for mask in range(1 << k):
        if sum(len(orbs[o]) for o in range(k) if mask >> o & 1) != 8:
            continue
        z1 = [o for o in range(k) if mask >> o & 1]
        z2 = [o for o in range(k) if not mask >> o & 1]
        for flips in itertools.product((0, 1), repeat=k):
            side = {}
            for o in z1:
                for u, bit in orbs[o].items():
                    side[u] = "AB"[bit ^ flips[o]]
            for o in z2:
                for u, bit in orbs[o].items():
                    side[u] = "CD"[bit ^ flips[o]]
            if all(side[ti[u]] == side[u] for u in range(N)):
                out.append(side)
    return out


def label_variants(side, ti):
    """All label assignments: within each side the two ti-twin pairs get X1/X2 in both orders."""
    pairs_by_side = {}
    for s in "ABCD":
        units = sorted(u for u in range(N) if side[u] == s)
        pairs = []
        seen = set()
        for u in units:
            if u in seen:
                continue
            pairs.append((u, ti[u]))
            seen |= {u, ti[u]}
        pairs_by_side[s] = pairs
    for bits in itertools.product((0, 1), repeat=4):
        lab = {}
        for s, bit in zip("ABCD", bits):
            p0, p1 = pairs_by_side[s]
            for u in p0:
                lab[u] = f"{s}{1 + bit}"
            for u in p1:
                lab[u] = f"{s}{2 - bit}"
        yield lab


def canonical_block_labelings(tj, tk, ti, L):
    """All label assignments of one block (tuple of 16 labels), reduced by the true label group
    acting on labels (units fixed): keep the lexicographic minimum of each orbit."""
    seen = set()
    for side in side_structures(tj, tk, ti):
        for lab in label_variants(side, ti):
            t = tuple(lab[u] for u in range(N))
            best = min(tuple(lam[POS[x]] for x in t) for lam in L)
            seen.add(best)
    return sorted(seen)


def incidence_certificate(design):
    """nauty certificate of the unit-port-hub incidence graph (colours: units / ports / hubs)."""
    # vertices: units 0..15, ports 16..63 (block b, l, r) -> 16 + 16b + 4l + r, hubs 64..87:
    # block b L-hub l -> 64 + 8b + l, R-hub r -> 64 + 8b + 4 + r
    adj = collections.defaultdict(set)
    for b in range(3):
        for l in range(4):
            for r in range(4):
                p = 16 + 16 * b + 4 * l + r
                adj[p] |= {64 + 8 * b + l, 64 + 8 * b + 4 + r}
                adj[64 + 8 * b + l].add(p)
                adj[64 + 8 * b + 4 + r].add(p)
    for u, t in enumerate(design):
        for b in range(3):
            for (l, r) in EDGES[t[b]]:
                p = 16 + 16 * b + 4 * l + r
                adj[u].add(p)
                adj[p].add(u)
    g = pynauty.Graph(88, adjacency_dict={k: sorted(v) for k, v in adj.items()},
                      vertex_coloring=[set(range(16)), set(range(16, 64)), set(range(64, 88))])
    return pynauty.certificate(g), pynauty.autgrp(g)[1]


def main():
    L = label_group()
    G1 = block1_stabiliser_true(L)
    print("label group:", len(L), "block-1 stabiliser:", len(G1), file=sys.stderr, flush=True)
    taus = valid_taus()
    seen = set()
    reps2 = []
    for t in taus:
        if t in seen:
            continue
        orbit = {conj(g, t) for g in G1}
        seen |= orbit
        reps2.append((t, len(orbit), [g for g in G1 if conj(g, t) == t]))
    print("tau_2 orbit representatives:", len(reps2), [r[1] for r in reps2], file=sys.stderr, flush=True)
    block1_labels = {}
    for s, units in SIDES1.items():
        block1_labels[units[0]] = block1_labels[units[1]] = f"{s}1"
        block1_labels[units[2]] = block1_labels[units[3]] = f"{s}2"
    classes = {}
    raw = 0
    for t2, orb2, stab2 in reps2:
        seen3 = set()
        for t3 in taus:
            if t3 in seen3:
                continue
            seen3 |= {conj(g, t3) for g in stab2}
            l2_list = canonical_block_labelings(TAU1, t3, t2, L)
            if not l2_list:
                continue
            l3_list = canonical_block_labelings(TAU1, t2, t3, L)
            if not l3_list:
                continue
            for lab2 in l2_list:
                for lab3 in l3_list:
                    design = [(block1_labels[u], lab2[u], lab3[u]) for u in range(N)]
                    raw += 1
                    cert, aut = incidence_certificate(design)
                    if cert not in classes:
                        classes[cert] = {"design": design, "aut": aut, "count": 0}
                    classes[cert]["count"] += 1
        print(f"tau_2 rep done: raw {raw}, classes {len(classes)}", file=sys.stderr, flush=True)
    out = []
    for cert, c in sorted(classes.items(), key=lambda kv: (-kv[1]["count"], kv[0])):
        d = c["design"]
        units = [{"labels": list(t), "ports": {str(i): [list(e) for e in EDGES[t[i]]] for i in range(3)}} for t in d]
        out.append({"raw_members": c["count"], "automorphisms": c["aut"], "design": d, "units": units})
    json.dump({"raw": raw, "classes": out}, sys.stdout, indent=1)
    print("raw designs:", raw, "isomorphism classes:", len(classes), file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
