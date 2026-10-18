#!/usr/bin/env python3
"""Self-contained construction of the (15,4) graph on 13608 vertices (route-chart lift with an odd-degree controller:
5-regular controller on 168 vertices, fibre F_3^4, alphabet a, a^-1, b, b^-1 and one self-inverse symbol s).
Writes the sorted edge list (vertices 0..13607, one "u v" per line, u < v) to the file given as argv[1].
No dependency on the search code.

Controller C: 168 vertices (c, g), c in {0,1,2,3}, g in Z_42, three permutations
  alpha(c, g) = (alpha0[c], g + va[c]),  alpha0 = (1, 2, 3, 0), va = (0, 9, 12, 0)
  beta (c, g) = (beta0[c],  g + vb[c]),  beta0  = (1, 2, 3, 0), vb = (19, 38, 7, 41)
  sigma(c, g) = (sigma0[c], g + vs[c]) for c < sigma0[c] and sigma = its own inverse on the other side,
                 sigma0 = (2, 3, 0, 1), vs[0] = 0 (link 0 <-> 2), vs[1] = 21 (link 1 <-> 3)
  i.e. sigma(0, g) = (2, g), sigma(2, g) = (0, g), sigma(1, g) = (3, g + 21), sigma(3, g) = (1, g - 21).
Directed edges of C carry the symbols a: v -> alpha(v), A = a^-1: v -> alpha^-1(v), b: v -> beta(v),
B = b^-1: v -> beta^-1(v), s: v -> sigma(v) (sigma is an involution: the s-edges form a perfect matching).
C is simple, 5-regular and has diameter 4.
Fibre F_3^4, x = (x0, x1, x2, x3) encoded as x0 + 3 x1 + 9 x2 + 27 x3; lift vertex (v, x) has index 81 v + x
with v = 42 c + g.
Chart: for each symbol t a matrix M_t in GL_4(F_3) (given by its columns M_t e_i) and a direction d_t:
  M_A = M_a^-1, d_A = M_a^-1 d_a; M_B = M_b^-1, d_B = M_b^-1 d_b (inverse pairs);
  M_s is an involution (M_s^2 = I) with M_s d_s = d_s (self-inverse symbol).
Edges: (v, x) ~ (t(v), M_t x + lambda d_t) for every symbol t in {a, A, b, B, s} and lambda in F_3.
Each vertex gets 5 x 3 = 15 neighbours, the relation is symmetric, order 168 x 81 = 13608, diameter 4.
"""
import sys

NG = 42
ALPHA0 = (1, 2, 3, 0); VA = (0, 9, 12, 0)
BETA0 = (1, 2, 3, 0);  VB = (19, 38, 7, 41)
SIGMA0 = (2, 3, 0, 1); VS = (0, 21, 0, 0)   # VS[c] used for c < SIGMA0[c]
CHART = {   # columns M_t e_0 .. e_3, direction d_t
    'a': ([[1, 0, 0, 0], [1, 1, 0, 1], [2, 0, 0, 2], [2, 1, 2, 0]], [0, 2, 0, 0]),
    'A': ([[1, 0, 0, 0], [0, 1, 1, 0], [2, 1, 1, 2], [2, 0, 2, 0]], [0, 2, 2, 0]),
    'b': ([[1, 1, 2, 0], [2, 1, 1, 2], [0, 0, 1, 1], [1, 0, 2, 1]], [2, 0, 0, 1]),
    'B': ([[2, 1, 1, 0], [0, 1, 0, 1], [1, 2, 1, 1], [2, 1, 0, 2]], [0, 0, 2, 2]),
    's': ([[1, 0, 0, 0], [1, 0, 1, 1], [0, 0, 1, 0], [2, 1, 2, 0]], [0, 1, 2, 1]),
}

def digits(x): return [(x // 3 ** i) % 3 for i in range(4)]
def undigits(d): return sum(d[i] * 3 ** i for i in range(4))
def matvec(cols, v): return [sum(v[i] * cols[i][j] for i in range(4)) % 3 for j in range(4)]

def main():
    alpha = {}; beta = {}; sigma = {}
    for c in range(4):
        for g in range(NG):
            alpha[(c, g)] = (ALPHA0[c], (g + VA[c]) % NG)
            beta[(c, g)] = (BETA0[c], (g + VB[c]) % NG)
            d = SIGMA0[c]
            sigma[(c, g)] = (d, (g + VS[c]) % NG) if c < d else (d, (g - VS[d]) % NG)
    ainv = {w: v for v, w in alpha.items()}; binv = {w: v for v, w in beta.items()}
    assert all(sigma[sigma[v]] == v for v in sigma)
    # chart consistency: inverse pairs and the involution
    for t, ti in (('a', 'A'), ('b', 'B')):
        for i in range(4):
            e = [1 if k == i else 0 for k in range(4)]
            assert matvec(CHART[ti][0], matvec(CHART[t][0], e)) == e
        assert matvec(CHART[ti][0], CHART[t][1]) == CHART[ti][1]
    for i in range(4):
        e = [1 if k == i else 0 for k in range(4)]
        assert matvec(CHART['s'][0], matvec(CHART['s'][0], e)) == e
    assert matvec(CHART['s'][0], CHART['s'][1]) == CHART['s'][1]
    vid = {(c, g): c * NG + g for c in range(4) for g in range(NG)}
    edges = set()
    for v in vid:
        for sym, tgt in (('a', alpha[v]), ('A', ainv[v]), ('b', beta[v]), ('B', binv[v]), ('s', sigma[v])):
            cols, d = CHART[sym]
            for x in range(81):
                Mx = matvec(cols, digits(x))
                for lam in range(3):
                    y = undigits([(Mx[j] + lam * d[j]) % 3 for j in range(4)])
                    u, w = vid[v] * 81 + x, vid[tgt] * 81 + y
                    assert u != w
                    edges.add((min(u, w), max(u, w)))
    with open(sys.argv[1], 'w') as f:
        for u, w in sorted(edges): f.write('%d %d\n' % (u, w))
    print('vertices', 4 * NG * 81, 'edges', len(edges))

main()
