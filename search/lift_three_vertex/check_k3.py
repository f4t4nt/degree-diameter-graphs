#!/usr/bin/env python3
"""Brute-force validation of gen_k3.py on a small group G(N, M, r, t) for a base spec:
  1. group axioms of the product (all triples), every map of A is an automorphism (all pairs);
  2. every symmetry used by the reduction is a lift isomorphism (explicit vertex bijection on random admissible
     assignments: automorphisms, loop inversion, permutations of parallel slots);
  3. orbit coverage: every orbit of admissible (simple) assignments under the group actually used by gen_k3 (A x loop
     inversions x parallel permutations) contains an emitted candidate; the over-coverage factor; also the orbit count
     under the full symmetry group including the shape automorphisms (for information).
usage: check_k3.py N M r t "SPEC" [shape-automorphisms as vertex permutations "012,021" (identity included or not)]
"""
import io
import itertools
import random
import sys

import numpy as np

import gen_k3


def main():
    N, M, r, t = (int(v) for v in sys.argv[1:5])
    spec = sys.argv[5]
    perms = [tuple(int(c) for c in p) for p in sys.argv[6].split(",")] if len(sys.argv) > 6 else []
    G = gen_k3.Group(N, M, r, t)
    m = G.m
    inv = G.inv.astype(np.int64)
    slots = gen_k3.parse_spec(spec)
    free = [s for s, sl in enumerate(slots) if not sl[3]]
    nfree = len(free)
    fails = 0
    # 1. axioms
    mult = np.zeros((m, m), dtype=np.int64)
    for a in range(m):
        for b in range(m):
            mult[a, b] = G.mul(a, b)
    ab_c = mult[mult, :]                                   # (a, b, c) -> (ab)c
    a_bc = mult[:, mult].transpose(0, 1, 2)                # (a, b, c) -> a(bc): mult[a, mult[b, c]]
    assoc_fail = int(np.sum(ab_c != a_bc))
    inv_fail = int(np.sum(mult[np.arange(m), inv] != 0) + np.sum(mult[inv, np.arange(m)] != 0))
    print(f"axioms: {assoc_fail} associativity failures on {m}^3 triples, {inv_fail} inverse failures")
    fails += assoc_fail + inv_fail
    nA = len(G.I0)
    img = np.zeros((nA, m), dtype=np.int64)
    for h in range(m):
        img[:, h] = G.images(G.I0, G.J0, G.L0, h)
    bad = 0
    for k in range(nA):
        p = img[k]
        if len(set(p.tolist())) != m or not np.array_equal(p[mult], mult[np.ix_(p, p)]):
            bad += 1
    print(f"automorphisms: |A| = {nA}, {bad} maps fail the homomorphism/bijection test")
    fails += bad

    # admissible values per free slot
    def allowed(s):
        kind = slots[s][0]
        if kind == "loop":
            return [x for x in range(1, m) if inv[x] != x]
        if kind == "semi":
            return [x for x in range(1, m) if inv[x] == x]
        if any(q[3] and q[:3] == slots[s][:3] for q in slots):
            return list(range(1, m))
        return list(range(m))

    ranges = [allowed(s) for s in free]
    total = 1
    for rg in ranges:
        total *= len(rg)
    if total == 0:
        print(f"RESULT N={N} M={M} r={r} t={t}: no admissible assignment (not enough involutions): nothing to check")
        return
    if total > 6e6:
        raise SystemExit(f"too many assignments to brute-force: {total}")
    # all assignments as a (total, nfree) array
    grids = np.meshgrid(*[np.array(rg, dtype=np.int64) for rg in ranges], indexing="ij")
    V = np.stack([gd.ravel() for gd in grids], axis=1)
    # simplicity: neighbours of (u, e) for each base vertex u
    simple = np.ones(len(V), dtype=bool)
    vals = {s: (V[:, free.index(s)] if s in free else np.zeros(len(V), dtype=np.int64)) for s in range(len(slots))}
    k = max(max(sl[1], sl[2]) for sl in slots) + 1
    for u in range(k):
        nb = []
        for s, (kind, a, b, fixed) in enumerate(slots):
            if kind == "link":
                if a == u:
                    nb.append(b * m + vals[s])
                if b == u:
                    nb.append(a * m + inv[vals[s]])
            elif kind == "loop":
                if a == u:
                    nb.append(u * m + vals[s])
                    nb.append(u * m + inv[vals[s]])
            elif a == u:
                nb.append(u * m + vals[s])
        NB = np.stack(nb, axis=1)
        simple &= ~np.any(NB == u * m, axis=1)
        NBs = np.sort(NB, axis=1)
        simple &= ~np.any(NBs[:, 1:] == NBs[:, :-1], axis=1)
    adm = V[simple]
    nadm = len(adm)
    radix = np.array([m] * nfree, dtype=np.int64)
    weights = np.cumprod(np.concatenate([[1], radix[:-1]]))[::-1].copy()
    weights = np.array([m ** (nfree - 1 - i) for i in range(nfree)], dtype=np.int64)
    codes = adm @ weights
    order = np.argsort(codes)
    scodes = codes[order]

    def index_of(W):
        c = W @ weights
        pos = np.searchsorted(scodes, c)
        pos = np.minimum(pos, len(scodes) - 1)
        ok = scodes[pos] == c
        return order[pos], ok

    # generators of the group used by gen_k3
    gens = []
    for kk in range(nA):
        gens.append(("aut", kk))
    loop_pos = [i for i, s in enumerate(free) if slots[s][0] == "loop"]
    for i in loop_pos:
        gens.append(("inv", i))
    groups = {}
    for i, s in enumerate(free):
        groups.setdefault(slots[s][:3], []).append(i)
    for pos in groups.values():
        for a in range(len(pos) - 1):
            gens.append(("swap", pos[a], pos[a + 1]))

    def apply(gen, W):
        W2 = W.copy()
        if gen[0] == "aut":
            W2 = img[gen[1]][W]
        elif gen[0] == "inv":
            W2[:, gen[1]] = inv[W[:, gen[1]]]
        else:
            W2[:, gen[1]], W2[:, gen[2]] = W[:, gen[2]], W[:, gen[1]]
        return W2

    # shape automorphisms (vertex permutations) as extra symmetries: the lift isomorphism (v, g) -> (perm v, g) moves
    # every slot value to the image slot (a link u->v whose image pair is stored as (v', u') with v' > u' is inverted,
    # the j-th slot of a kind at a place goes to the j-th slot at the image place), after which the spanning tree may
    # carry non-identity voltages: re-gauge with h_0 = e, h_v = w^-1 h_u along the tree dart u->v (voltage w), every
    # dart u->v getting h_u^-1 w h_v (loops / semi-edges at u are conjugated by h_u); the lift isomorphism is then
    # (v, g) -> (perm v, g h_v)
    nslots = len(slots)
    place_all = {}
    for s in range(nslots):
        place_all.setdefault(slots[s][:3], []).append(s)
    tree = [s for s in range(nslots) if slots[s][3]]

    def apply_perm(perm, W, return_h=False):
        full = np.zeros((len(W), nslots), dtype=np.int64)
        for i, s in enumerate(free):
            full[:, s] = W[:, i]
        img2 = np.zeros_like(full)
        used = {}
        for s in range(nslots):
            kind, a, b, fixed = slots[s]
            if kind == "link":
                a2, b2 = perm[a], perm[b]
                invert = a2 > b2
                key = ("link", min(a2, b2), max(a2, b2))
            else:
                invert = False
                key = (kind, perm[a], perm[a])
            if key not in place_all or used.get(key, 0) >= len(place_all[key]):
                raise SystemExit(f"permutation {perm} is not a shape automorphism (slot {s} -> {key})")
            j = place_all[key][used.get(key, 0)]
            used[key] = used.get(key, 0) + 1
            img2[:, j] = inv[full[:, s]] if invert else full[:, s]
        # re-gauge along the tree (rooted at vertex 0; tree darts are oriented u < v in the spec)
        h = {0: np.zeros(len(W), dtype=np.int64)}
        pending = list(tree)
        while pending:
            progress = False
            for s in list(pending):
                kind, u, v, fixed = slots[s]
                if u in h and v not in h:
                    h[v] = mult[inv[img2[:, s]], h[u]]
                elif v in h and u not in h:
                    h[u] = mult[img2[:, s], h[v]]
                else:
                    continue
                pending.remove(s)
                progress = True
            if not progress:
                raise SystemExit("tree is not connected")
        out = np.zeros_like(img2)
        for s in range(nslots):
            kind, u, v, fixed = slots[s]
            if kind == "link":
                out[:, s] = mult[mult[inv[h[u]], img2[:, s]], h[v]]
            else:
                out[:, s] = mult[mult[inv[h[u]], img2[:, s]], h[u]]
        if np.any(out[:, tree] != 0):
            raise SystemExit("re-gauging failed")
        return (out[:, free], h) if return_h else out[:, free]

    # 2. symmetries as lift isomorphisms on random admissible assignments
    rng = random.Random(7)

    def lift_edges(v):
        vals1 = {s: 0 for s in range(len(slots))}
        for s, x in zip(free, v):
            vals1[s] = int(x)
        E = set()
        for gg in range(m):
            for s, (kind, u, w, fixed) in enumerate(slots):
                h = int(mult[gg, vals1[s]])
                E.add(frozenset(((u, gg), (w, h))) if kind == "link" else frozenset(((u, gg), (u, h))))
        return E

    nsym = 0
    if nadm == 0:
        print(f"RESULT N={N} M={M} r={r} t={t}: no admissible (simple) assignment exists: nothing to check")
        return
    for _ in range(6):
        v = adm[rng.randrange(nadm)]
        E1 = lift_edges(v)
        tests = rng.sample(gens, min(5, len(gens)))
        for gen in tests:
            w = apply(gen, v[None, :])[0]
            if gen[0] == "aut":
                p = img[gen[1]]
                psi = lambda x, p=p: (x[0], int(p[x[1]]))
            else:
                psi = lambda x: x
            nsym += 1
            if {frozenset(psi(x) for x in e) for e in E1} != lift_edges(w):
                fails += 1
                print(f"  symmetry failed: {gen} {v} -> {w}")
        for perm in perms:
            w, hh = apply_perm(perm, v[None, :], return_h=True)
            w = w[0]
            nsym += 1
            if {frozenset((perm[x[0]], int(mult[x[1], hh[perm[x[0]]][0]])) for x in e) for e in E1} != lift_edges(w):
                fails += 1
                print(f"  shape automorphism failed: {perm} {v} -> {w}")
    print(f"lift isomorphisms: {nsym} checks, failures so far {fails}")

    # 3. orbit coverage by label propagation
    def orbits(gen_list, extra_perms):
        lab = np.arange(nadm)
        maps = []
        for gen in gen_list:
            idx, ok = index_of(apply(gen, adm))
            if not np.all(ok):
                raise SystemExit(f"generator {gen} maps an admissible assignment outside the admissible set")
            maps.append(idx)
        for perm in extra_perms:
            idx, ok = index_of(apply_perm(perm, adm))
            if not np.all(ok):
                raise SystemExit(f"shape automorphism {perm} maps outside the admissible set")
            maps.append(idx)
        while True:
            new = lab.copy()
            for idx in maps:
                new = np.minimum(new, new[idx])
                np.minimum.at(new, idx, new)
            if np.array_equal(new, lab):
                break
            lab = new
        return lab

    lab = orbits(gens, [])
    buf = io.StringIO()
    lines, cands = gen_k3.generate(G, spec, False, buf)
    emitted = []
    for ln in buf.getvalue().splitlines():
        p = [int(q) for q in ln.split()]
        if len(p) >= nfree and p[nfree - 1] == -1:          # compact line: prefix, -1, explicit last-slot values
            for last in p[nfree:]:
                emitted.append(p[:nfree - 1] + [last])
        elif len(p) == nfree:
            emitted.append(p)
        else:
            for last in range(p[-2], p[-1]):
                emitted.append(p[:-2] + [last])
    E = np.array(emitted, dtype=np.int64).reshape(-1, nfree)
    idx, ok = index_of(E)
    eidx = idx[ok]
    covered = set(lab[eidx].tolist())
    all_orbits = set(np.unique(lab).tolist())
    missing = sorted(all_orbits - covered)
    print(f"orbit coverage: {len(all_orbits)} orbits of {nadm} admissible assignments, {len(eidx)} emitted admissible "
          f"candidates ({len(E)} emitted, {cands} counted), {len(missing)} orbits MISSED, over-coverage "
          f"{len(eidx) / max(1, len(all_orbits)):.2f}")
    if missing:
        fails += len(missing)
        print("  missed e.g.", [adm[np.nonzero(lab == q)[0][0]].tolist() for q in missing[:5]])
    if perms:
        lab2 = orbits(gens, perms)
        print(f"orbits under the full symmetry group (with the shape automorphisms {perms}): {len(np.unique(lab2))}")
    print(f"RESULT N={N} M={M} r={r} t={t} spec={spec}: {'OK' if fails == 0 else 'FAIL'} ({fails} failures)")


if __name__ == "__main__":
    main()
