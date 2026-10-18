#!/usr/bin/env python3
"""Self-contained construction of the (15,3) graph on 1620 vertices (route-chart lift with 27-point cells over a
5-regular 60-vertex controller).  Writes the sorted edge list (vertices 0..1619, one "u v" per line, u < v) to the
file given as argv[1].  No dependency on the search code; the per-dart chart is read from charts_1620_15_3.json
beside this script.

Controller: the regular Z_20-lift of the 3-vertex Schreier base alpha0 = (0,2,1) [a fixes base vertex 0, swaps 1
  and 2], beta0 = (1,2,0) [3-cycle], sigma0 = (1,0,2) [link 0-1, semi-edge at 2] with voltages va = 6,3,12
  (a-darts: (c, g) -> (alpha0 c, g + va[c])), vb = 0,0,0, vs = 1 (link (0, g) -> (1, g + 1)), 10 (semi-edge
  (2, g) -> (2, g + 10)); vertex index c * 20 + g.  Every vertex is on a triangle and every non-adjacent ordered
  pair is joined by a walk of length 3 (A^3 + A > 0), the diameter-3 controller condition.
Cells F_3^3 (27 points), vertex index 27 v + 9 x0 + 3 x1 + x2.  Each dart v -> w (slot i of v) carries a matrix M in
  GL(3,3) and a direction d (the reverse dart carries M^-1 and M^-1 d), and (v, x) ~ (w, M x + lambda d) for
  lambda in F_3.  The chart is of the coplanar type (four departure directions of each vertex on a projective line
  of PG(2,3), one special); the per-edge types come from the type SAT (typesat3.py) and the matrices from its
  brute-force realization.
"""
import itertools
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CHARTS = os.path.join(HERE, 'charts_1620_15_3.json')


def controller():
    """the 60-vertex controller from its Schreier description (independent of the chart file's adjacency)."""
    nG, k = 20, 3
    alpha0, beta0, sigma0 = (0, 2, 1), (1, 2, 0), (1, 0, 2)
    va, vb = (6, 3, 12), (0, 0, 0)
    vs_link, vs_semi = 1, 10        # link 0 <-> 1 voltage 1; semi-edge at base vertex 2 voltage 10
    N = k * nG
    adj = [[None] * 5 for _ in range(N)]
    alpha = {}; beta = {}; sig = {}
    for c in range(k):
        for g in range(nG):
            v = c * nG + g
            alpha[v] = alpha0[c] * nG + (g + va[c]) % nG
            beta[v] = beta0[c] * nG + (g + vb[c]) % nG
            if c == 0: sig[v] = 1 * nG + (g + vs_link) % nG
            elif c == 1: sig[v] = 0 * nG + (g - vs_link) % nG
            else: sig[v] = 2 * nG + (g + vs_semi) % nG
    ainv = {alpha[v]: v for v in range(N)}; binv = {beta[v]: v for v in range(N)}
    for v in range(N):
        adj[v] = [alpha[v], ainv[v], beta[v], binv[v], sig[v]]
    return N, adj


def main():
    chart = json.load(open(CHARTS))
    N, adj = controller()
    assert N == chart['N'] and adj == chart['adj'], 'controller mismatch with the chart file'
    pts = list(itertools.product(range(3), repeat=3)); idx = {p: i for i, p in enumerate(pts)}
    def mv(M, x): return tuple(sum(M[i][j] * x[j] for j in range(3)) % 3 for i in range(3))
    edges = set()
    for v in range(N):
        for i in range(5):
            dart = chart['darts'][f'{v},{i}']
            w, M, d = dart['to'], dart['M'], dart['d']
            assert w == adj[v][i]
            for x in pts:
                y0 = mv(M, x)
                for lam in range(3):
                    y = tuple((y0[t] + lam * d[t]) % 3 for t in range(3))
                    a, b = v * 27 + idx[x], w * 27 + idx[y]
                    assert a != b
                    edges.add((min(a, b), max(a, b)))
    with open(sys.argv[1], 'w') as f:
        for a, b in sorted(edges): f.write(f'{a} {b}\n')
    print('vertices', 27 * N, 'edges', len(edges))


main()
