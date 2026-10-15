"""Brute-force validation of gen_b2.py on a small group G(N, M, r, t), for any shape in gen_b2.SPECS:
  1. group axioms of the bitlift product (associativity on all triples, inverses), and every map of A is an
     automorphism (checked on all pairs);
  2. every symmetry used by the reduction is a lift isomorphism (explicit vertex bijection on random assignments:
     automorphisms, loop inversion, permutations of parallel slots, the B2a vertex swap), plus the unused extras
     (B2a link swap, Haar vertex swap / re-gauging) for information;
  3. every orbit of admissible (simple) assignments under the symmetry group actually used by gen_b2 contains an
     emitted candidate (orbit coverage), and the over-coverage factor.
usage: check_b2.py N M r t SHAPE
"""
import io
import itertools
import random
import sys

import numpy as np

import gen_b2


def lift_edges(G, shape, v):
    """edge set of the lift of SPECS[shape] with the free-slot values v (fixed links carry e); semi-edge: (u,g)~(u,g s)"""
    m = G.m
    E = set()
    slots = gen_b2.parse_spec(gen_b2.SPECS[shape])
    free = [s for s, sl in enumerate(slots) if not sl[3]]
    vals = {s: 0 for s in range(len(slots))}
    for s, x in zip(free, v):
        vals[s] = x
    for g in range(m):
        for s, (kind, u, w, fixed) in enumerate(slots):
            h = G.mul(g, vals[s])
            if kind == "link":
                E.add(frozenset(((u, g), (w, h))))
            else:
                E.add(frozenset(((u, g), (u, h))))
    return E


def main():
    N, M, r, t = (int(v) for v in sys.argv[1:5])
    shape = sys.argv[5]
    G = gen_b2.Group(N, M, r, t)
    m = G.m
    inv = [int(v) for v in G.inv]
    slots = gen_b2.parse_spec(gen_b2.SPECS[shape])
    free = [s for s, sl in enumerate(slots) if not sl[3]]
    fails = 0
    # 1. axioms
    for a in range(m):
        for b in range(m):
            ab = G.mul(a, b)
            for c in range(m):
                if G.mul(ab, c) != G.mul(a, G.mul(b, c)):
                    fails += 1
    for a in range(m):
        if G.mul(a, inv[a]) != 0 or G.mul(inv[a], a) != 0 or G.mul(a, 0) != a or G.mul(0, a) != a:
            fails += 1
    print(f"axioms: {fails} failures on {m}^3 triples")
    img = np.zeros((len(G.I0), m), dtype=np.int64)
    for h in range(m):
        img[:, h] = G.images(G.I0, G.J0, G.L0, h)
    mult = np.zeros((m, m), dtype=np.int64)
    for a in range(m):
        for b in range(m):
            mult[a, b] = G.mul(a, b)
    bad = 0
    for k in range(len(G.I0)):
        p = img[k]
        if len(set(p.tolist())) != m or not np.array_equal(p[mult], mult[np.ix_(p, p)]):
            bad += 1
    print(f"automorphisms: |A| = {len(G.I0)}, {bad} maps fail the homomorphism/bijection test")
    fails += bad
    # admissible values per free slot and simplicity (distinct neighbours of (u, e) for each u)
    loops = [g for g in range(1, m) if inv[g] != g]
    invols = [g for g in range(1, m) if inv[g] == g]

    def allowed(s):
        kind = slots[s][0]
        if kind == "loop":
            return loops
        if kind == "semi":
            return invols
        if any(q[3] and q[:3] == slots[s][:3] for q in slots):
            return list(range(1, m))
        return list(range(m))

    def simple(v):
        vals = {s: 0 for s in range(len(slots))}
        for s, x in zip(free, v):
            vals[s] = x
        for u in (0, 1):
            nb = []
            for s, (kind, a, b, fixed) in enumerate(slots):
                if kind == "link":
                    if a == u:
                        nb.append((b, vals[s]))
                    if b == u:
                        nb.append((a, inv[vals[s]]))
                elif kind == "loop":
                    if a == u:
                        nb += [(u, vals[s]), (u, inv[vals[s]])]
                elif a == u:
                    nb.append((u, vals[s]))
            if (u, 0) in nb or len(set(nb)) != len(nb):
                return False
        return True

    # parallel groups of free slots (same kind, same place) -> permutations are symmetries
    groups = {}
    for i, s in enumerate(free):
        groups.setdefault(slots[s][:3], []).append(i)
    loop_pos = [i for i, s in enumerate(free) if slots[s][0] == "loop"]

    # vertex swap (symmetric shapes, gen_b2.symmetric_rule): the j-th free slot of kind K at vertex u takes the value of
    # the j-th free slot of kind K at vertex 1 - u, links are inverted (the tree link e stays e)
    sym = bool(gen_b2.symmetric_rule(slots))

    def swapped(v):
        w = list(v)
        for i, s in enumerate(free):
            kind, u, _, _ = slots[s]
            if kind == "link":
                w[i] = inv[v[i]]
            else:
                same = [q for q in range(len(free)) if slots[free[q]][0] == kind and slots[free[q]][1] == u]
                other = [q for q in range(len(free)) if slots[free[q]][0] == kind and slots[free[q]][1] == 1 - u]
                w[i] = v[other[same.index(i)]]
        return tuple(w)

    def gens_used(v):
        out = [tuple(int(img[k][x]) for x in v) for k in range(len(G.I0))]
        for i in loop_pos:
            w = list(v); w[i] = inv[v[i]]; out.append(tuple(w))
        for pos in groups.values():
            if len(pos) > 1:
                for perm in itertools.permutations(pos):
                    if perm != tuple(pos):
                        w = list(v)
                        for a, b in zip(pos, perm):
                            w[a] = v[b]
                        out.append(tuple(w))
        if sym:
            out.append(swapped(v))
        return out

    def gens_extra(v):
        if shape == "B2a":
            a, b, c = v
            return [(a, G.mul(G.mul(c, b), inv[c]), inv[c])]
        if shape == "haar":
            c1 = v[0]
            return [tuple(inv[x] for x in v), (inv[c1], G.mul(v[1], inv[c1]), G.mul(v[2], inv[c1]))]
        return []

    # 2. symmetries as lift isomorphisms on random admissible assignments
    rng = random.Random(7)

    def psi_for(v, w):
        """vertex bijection realising v -> w for the generators above (identity for inversions/permutations)"""
        for k in range(len(G.I0)):
            if tuple(int(img[k][x]) for x in v) == w:
                p = img[k]
                return lambda x: (x[0], int(p[x[1]]))
        if sym and w == swapped(v):
            return lambda x: (1 - x[0], x[1])
        if shape == "B2a" and w == (v[0], G.mul(G.mul(v[2], v[1]), inv[v[2]]), inv[v[2]]):
            c = v[2]
            return lambda x: x if x[0] == 0 else (1, G.mul(x[1], inv[c]))
        if shape == "haar" and w == tuple(inv[x] for x in v):
            return lambda x: (1 - x[0], x[1])
        if shape == "haar" and w == (inv[v[0]], G.mul(v[1], inv[v[0]]), G.mul(v[2], inv[v[0]])):
            c1 = v[0]
            return lambda x: x if x[0] == 0 else (1, G.mul(x[1], inv[c1]))
        return lambda x: x

    nsym = 0
    for _ in range(8 if all(allowed(s) for s in free) else 0):
        v = None
        for _try in range(2000):   # shapes needing more involutions than the group has may have no simple assignment
            w = tuple(rng.choice(allowed(s)) for s in free)
            if simple(w):
                v = w
                break
        if v is None:
            break
        E1 = lift_edges(G, shape, v)
        tests = rng.sample(gens_used(v), min(6, len(gens_used(v)))) + gens_extra(v)
        for w in tests:
            psi = psi_for(v, w)
            nsym += 1
            if {frozenset(psi(x) for x in e) for e in E1} != lift_edges(G, shape, w):
                fails += 1
                print(f"  symmetry failed: {v} -> {w}")
    print(f"lift isomorphisms: {nsym} checks, failures so far {fails}")
    # 3. orbit coverage
    buf = io.StringIO()
    lines, cands = gen_b2.generate(G, shape, False, buf)
    emitted = set()
    for ln in buf.getvalue().splitlines():
        p = [int(q) for q in ln.split()]
        if len(p) == len(free):          # full line (orbit representative of the last slot)
            emitted.add(tuple(p))
        else:                            # "v1 .. v_{d} lo hi": range for the last slot
            for last in range(p[-2], p[-1]):
                emitted.add(tuple(p[:-2]) + (last,))
    adm = [v for v in itertools.product(*(allowed(s) for s in free)) if simple(v)]
    index = {v: i for i, v in enumerate(adm)}
    parent = list(range(len(adm)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        x, y = find(x), find(y)
        if x != y:
            parent[x] = y

    for v in adm:
        for w in gens_used(v):
            union(index[v], index[w])
    roots = {}
    for v in adm:
        roots.setdefault(find(index[v]), []).append(v)
    emitted_adm = [v for v in emitted if v in index]
    covered = {find(index[v]) for v in emitted_adm}
    missing = [orb[0] for root, orb in roots.items() if root not in covered]
    print(f"orbit coverage ({shape}, group used by gen_b2): {len(roots)} orbits of {len(adm)} admissible assignments, "
          f"{len(emitted_adm)} emitted admissible candidates ({len(emitted)} emitted, {cands} counted), "
          f"{len(missing)} orbits MISSED, over-coverage {len(emitted_adm) / max(1, len(roots)):.2f}")
    if missing:
        fails += len(missing)
        print("  missed e.g.", missing[:5])
    for v in adm:
        for w in gens_extra(v):
            if w in index:
                union(index[v], index[w])
    full = len({find(index[v]) for v in adm})
    print(f"orbits under the full symmetry group (with the unused extra symmetries): {full}")
    print(f"RESULT N={N} M={M} r={r} t={t} shape={shape}: {'OK' if fails == 0 else 'FAIL'} ({fails} failures)")


if __name__ == "__main__":
    main()
