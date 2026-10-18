#!/usr/bin/env python3
"""Self-contained construction of the (12,4) graph on 6156 vertices (route-chart lift, Mizuno's fibre chart over
F_3^4 on a 76-vertex Z_19-lift controller).  Writes the edge list (vertices 0..6155) to the file given as argv[1].

Controller C (76 vertices (c, g), c in {0,1,2,3}, g in Z_19), two permutations:
  alpha(c, g) = (c+1 mod 4, g + va[c]),  va = (0, 12, 4, 0)
  beta (c, g) = (beta0[c],   g + vb[c]),  beta0 = (1, 3, 0, 2), vb = (17, 0, 0, 1)
Directed edges of C carry the symbols a: v -> alpha(v), A: v -> alpha^-1(v), b: v -> beta(v), B: v -> beta^-1(v).
Fibre F = F_3^4, x = (x0, x1, x2, x3) encoded as x0 + 3 x1 + 9 x2 + 27 x3.  Each symbol t has a linear map
M_t (columns M_t e_i below, defined modulo the translation line T_t) and a translation line T_t = F_3 * d_t.
Lift vertex (v, x), v = 4*... index = v*81 + x, is adjacent to (t(v), M_t x + lambda d_t) for every symbol t and
lambda in F_3.  Degree 4 * 3 = 12, order 76 * 81 = 6156, diameter 4.
"""
import sys

VA = (0, 12, 4, 0); VB = (17, 0, 0, 1); BETA0 = (1, 3, 0, 2); NG = 19
# chart (columns M_t e_0..e_3, direction d_t), read off R. Mizuno's published 5832-vertex (12,4) graph
CHART = {
    'a': ([[0, 0, 1, 1], [0, 1, 1, 1], [0, 1, 1, 0], [0, 0, 0, 0]], [1, 0, 0, 0]),
    'A': ([[0, 0, 0, 0], [2, 1, 0, 0], [1, 2, 1, 0], [0, 1, 2, 0]], [0, 0, 0, 1]),
    'b': ([[2, 0, 1, 2], [1, 0, 0, 1], [0, 0, 0, 0], [0, 0, 1, 2]], [0, 1, 0, 0]),
    'B': ([[2, 0, 0, 1], [0, 0, 0, 0], [1, 1, 0, 0], [1, 1, 0, 2]], [0, 0, 1, 0]),
}

def digits(x): return [(x // 3 ** i) % 3 for i in range(4)]
def undigits(d): return sum(d[i] * 3 ** i for i in range(4))

def main():
    N = 4 * NG
    alpha = {}; beta = {}
    for c in range(4):
        for g in range(NG):
            alpha[(c, g)] = ((c + 1) % 4, (g + VA[c]) % NG)
            beta[(c, g)] = (BETA0[c], (g + VB[c]) % NG)
    ainv = {w: v for v, w in alpha.items()}; binv = {w: v for v, w in beta.items()}
    vid = {(c, g): c * NG + g for c in range(4) for g in range(NG)}
    edges = set()
    for v in vid:
        for sym, tgt in (('a', alpha[v]), ('A', ainv[v]), ('b', beta[v]), ('B', binv[v])):
            cols, d = CHART[sym]
            for x in range(81):
                xd = digits(x)
                Mx = [sum(xd[i] * cols[i][j] for i in range(4)) % 3 for j in range(4)]
                for lam in range(3):
                    y = undigits([(Mx[j] + lam * d[j]) % 3 for j in range(4)])
                    u, w = vid[v] * 81 + x, vid[tgt] * 81 + y
                    edges.add((min(u, w), max(u, w)))
    with open(sys.argv[1], 'w') as f:
        for u, w in sorted(edges): f.write('%d %d\n' % (u, w))
    print('vertices', N * 81, 'edges', len(edges))

main()
