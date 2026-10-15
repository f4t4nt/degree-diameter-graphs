#!/usr/bin/env python3
"""make_edges.py: independent rebuild (gf.py's group law from the file's JSON spec, no use of table_sweep.c's table) of the edge list
of a candidate found by table_sweep, written sorted as "u v" (u < v) lines.

  make_edges.py GRP cay7 a b c OUT.edges        Cay(G, {y^(M/2)} u {a, b, c}^+-), vertex = element index k*M + m
  make_edges.py GRP lift a b c2 c3 OUT.edges    shape-103 lift: vertex (v, g) -> v*|G| + g; (0,g)~(0,g a^+-1), (1,g)~(1,g b^+-1),
                                                (0,g)~(1,g), (1,g c2), (1,g c3)
"""
import sys

import numpy as np

import gf  # noqa: E402


def main():
    gfile, fam = sys.argv[1], sys.argv[2]
    nK, M, TK, alpha, spec = gf.read_grp(gfile)
    blocks = [gf.build_block(bl) for bl in spec["blocks"]]
    nK2, TK2, alpha2, gens = gf.direct_product(blocks)
    assert nK2 == nK and np.array_equal(TK2, TK) and np.array_equal(alpha2, alpha), "group file does not match its spec"
    G = gf.Group(nK, M, TK2, alpha2, spec)
    n = G.n
    ALL = np.arange(n)
    inv = G.inverse_all()
    edges = set()
    if fam == "cay7":
        a, b, c = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]); out = sys.argv[6]
        t = M // 2
        assert inv[t] == t and t != 0
        for s in [t, a, inv[a], b, inv[b], c, inv[c]]:
            W = G.mul(ALL, np.full(n, s))
            for u, w in zip(ALL.tolist(), W.tolist()):
                if u != w:
                    edges.add((min(u, w), max(u, w)))
    elif fam == "lift":
        a, b, c2, c3 = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]); out = sys.argv[7]
        for s in [a, inv[a]]:
            W = G.mul(ALL, np.full(n, s))
            for u, w in zip(ALL.tolist(), W.tolist()):
                edges.add((min(u, w), max(u, w)))
        for s in [b, inv[b]]:
            W = G.mul(ALL, np.full(n, s))
            for u, w in zip(ALL.tolist(), W.tolist()):
                edges.add((n + min(u, w), n + max(u, w)))
        for s in [0, c2, c3]:
            W = G.mul(ALL, np.full(n, s))
            for u, w in zip(ALL.tolist(), W.tolist()):
                edges.add((u, n + w))
    else:
        raise SystemExit("family cay7|lift")
    with open(out, "w") as f:
        for u, v in sorted(edges):
            f.write(f"{u} {v}\n")
    print(f"edges: {len(edges)} vertices: {n if fam == 'cay7' else 2 * n} -> {out}")


if __name__ == "__main__":
    main()
