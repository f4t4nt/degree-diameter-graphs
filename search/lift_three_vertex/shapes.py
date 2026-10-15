#!/usr/bin/env python3
"""Connected 3-vertex degree-5 base pregraphs (loops l_v, semi-edges s_v, link multiplicities m_uv with
2 l_v + s_v + sum_u m_uv = 5) up to isomorphism, with their automorphism group order, the free-slot structure and a
bitlift base spec (k=3;deg=5;...; `*` = spanning-tree link fixed to the identity).

Slot order in the spec (the sequential-stabiliser generator reduces the first slots the most and leaves the last free
slot as a range): semi-edges first (involutions: the smallest ranges), then loops, then the free (non-tree) links;
within a kind, vertex 0, 1, 2; parallel slots adjacent.  The tree is the two link pairs (0,1), (0,2) when both exist,
else the two existing pairs.

usage: shapes.py            prints the table (one line per shape) and the count
       shapes.py --json     one JSON object per line: id, m01, m02, m12, loops, semis, aut, spec, free slots, S
"""
import itertools
import json
import sys


def enumerate_shapes():
    seen = {}
    pairs = [(0, 1), (0, 2), (1, 2)]
    for m01 in range(6):
        for m02 in range(6):
            for m12 in range(6):
                mult = {(0, 1): m01, (0, 2): m02, (1, 2): m12}
                # connected on 3 vertices: at least two of the three pairs present
                if sum(1 for p in pairs if mult[p] > 0) < 2:
                    continue
                ldeg = [m01 + m02, m01 + m12, m02 + m12]
                if max(ldeg) > 5:
                    continue
                opts = []
                for v in range(3):
                    rem = 5 - ldeg[v]
                    opts.append([(l, rem - 2 * l) for l in range(rem // 2 + 1)])
                for combo in itertools.product(*opts):
                    loops = tuple(c[0] for c in combo)
                    semis = tuple(c[1] for c in combo)
                    key, aut = canon(mult, loops, semis)
                    if key not in seen:
                        seen[key] = (mult, loops, semis, aut)
    return list(seen.values())


def canon(mult, loops, semis):
    """canonical key over the 6 vertex permutations; aut = number of permutations fixing the labelled shape"""
    best = None
    aut = 0
    me = None
    for perm in itertools.permutations(range(3)):
        pm = {}
        for (u, v), m in mult.items():
            a, b = sorted((perm[u], perm[v]))
            pm[(a, b)] = m
        pl = [0] * 3
        ps = [0] * 3
        for v in range(3):
            pl[perm[v]] = loops[v]
            ps[perm[v]] = semis[v]
        key = (pm[(0, 1)], pm[(0, 2)], pm[(1, 2)], tuple(pl), tuple(ps))
        if me is None:
            me = key
        if key == me:
            aut += 1
        if best is None or key < best:
            best = key
    return best, aut


def spec_of(mult, loops, semis):
    """bitlift spec string + the list of free slots (kind, u, v)"""
    tree = [p for p in ((0, 1), (0, 2), (1, 2)) if mult[p] > 0][:2]
    toks = ["k=3", "deg=5"]
    free = []
    for v in range(3):
        for _ in range(semis[v]):
            toks.append(f"semi:{v}")
            free.append(("semi", v, v))
    for v in range(3):
        for _ in range(loops[v]):
            toks.append(f"loop:{v}")
            free.append(("loop", v, v))
    for p in ((0, 1), (0, 2), (1, 2)):
        m = mult[p]
        if m == 0:
            continue
        if p in tree:
            toks.append(f"link*:{p[0]},{p[1]}")
            m -= 1
        for _ in range(m):
            toks.append(f"link:{p[0]},{p[1]}")
            free.append(("link", p[0], p[1]))
    return ";".join(toks), free


def main():
    shapes = enumerate_shapes()
    rows = []
    for i, (mult, loops, semis, aut) in enumerate(sorted(shapes, key=lambda s: (sum(s[2]), sum(s[1]), -s[0][(0, 1)], -s[0][(0, 2)], -s[0][(1, 2)]))):
        spec, free = spec_of(mult, loops, semis)
        rows.append({"id": f"s{i:02d}", "m01": mult[(0, 1)], "m02": mult[(0, 2)], "m12": mult[(1, 2)], "loops": list(loops),
                     "semis": list(semis), "S": sum(semis), "P": sum(loops), "L": sum(mult.values()), "aut": aut,
                     "nfree": len(free), "free": free, "spec": spec})
    if "--json" in sys.argv:
        for r in rows:
            print(json.dumps(r))
        return
    print("id   m01 m02 m12  loops   semis   S P L aut nfree  spec")
    for r in rows:
        print(f"{r['id']}  {r['m01']}   {r['m02']}   {r['m12']}   {r['loops']} {r['semis']} {r['S']} {r['P']} {r['L']} {r['aut']:3d} {r['nfree']:4d}   {r['spec']}")
    print(f"# {len(rows)} connected 3-vertex degree-5 pregraphs up to isomorphism; semi-edge totals: "
          + ", ".join(f"S={s}: {sum(1 for r in rows if r['S'] == s)}" for s in sorted({r['S'] for r in rows})))


if __name__ == "__main__":
    main()
