#!/usr/bin/env python3
"""Self-contained construction of the (14,4) graph on 9720 vertices: a route-chart lift with the six-symbol MIXED
alphabet over F_3^4 (four F_3-line symbols a, a^-1, b, b^-1 with 3 neighbours each, one plain inverse pair p, p^-1
with 1 neighbour each) of the Cayley graph of G = SL(2,3) x Z_5 (order 120).  No dependency on the search code.
Writes the sorted edge list (vertices 0..9719, one "u v" per line, u < v) to the file given as argv[1].

Group G = SL(2,3) x Z_5: pairs (X, i), X a 2x2 matrix over F_3 of determinant 1 (24 of them), i in Z_5,
(X, i) (Y, j) = (X Y mod 3, i + j mod 5).  Element numbering: SL(2,3) listed as the identity first, then the
remaining matrices ((a, b), (c, d)) with det 1 in lexicographic order of (a, b, c, d); G element (X, i) gets the
index 5 * index(X) + i (so the identity is 0 and p = (I, 1) is 1).
Generators (indices in that numbering): a = 16 = ([[0, 1], [2, 2]], 1) of order 15, b = 12 = ([[0, 1], [2, 1]], 2)
of order 30, plain pair p = 1 = (I, 1) (central, order 5), p^-1 = 4.
Controller C = Cay(G, {a, a^-1, b, b^-1, p, p^-1}): g -> g a, g a^-1, g b, g b^-1, g p, g p^-1 (right multiplication),
a 6-regular graph on 120 vertices; the four line symbols form a 4-regular graph of diameter 7 on their own, p is
central so the plain edges commute with the line edges.
Fibre F_3^4, x = (x0, x1, x2, x3) encoded as x0 + 3 x1 + 9 x2 + 27 x3; lift vertex (g, x) has index 81 g + x.
Chart: for each line symbol t a matrix A_t in GL_4(F_3) (rows below, y_i = sum_j A[i][j] x_j) and a direction d_t
with A_{a^-1} = A_a^-1, d_{a^-1} = A_a^-1 d_a (same for b); for the plain pair M_p = -I, c_p = (2, 0, 2, 2), so that
M_{p^-1} = M_p^-1 = -I and c_{p^-1} = -M_p^-1 c_p = c_p.
Edges: (g, x) ~ (g t, A_t x + lambda d_t) for t in {a, a^-1, b, b^-1}, lambda in F_3 (12 neighbours), and
(g, x) ~ (g p, -x + c_p), (g, x) ~ (g p^-1, -x + c_p) (2 neighbours): 14-regular, 9720 vertices, diameter 4.
"""
import sys

# --- the group ---------------------------------------------------------------------------------------------------
def mat_mul(X, Y):
    (a, b), (c, d) = X; (e, f), (g, h) = Y
    return (((a * e + b * g) % 3, (a * f + b * h) % 3), ((c * e + d * g) % 3, (c * f + d * h) % 3))

I2 = ((1, 0), (0, 1))
SL = [((a, b), (c, d)) for a in range(3) for b in range(3) for c in range(3) for d in range(3) if (a * d - b * c) % 3 == 1]
SL.remove(I2); SL.insert(0, I2)
ELEMS = [(X, i) for X in SL for i in range(5)]
INDEX = {g: k for k, g in enumerate(ELEMS)}
def gmul(u, v): return (mat_mul(u[0], v[0]), (u[1] + v[1]) % 5)
def ginv(u):
    for v in ELEMS:
        if gmul(u, v) == (I2, 0): return v
    raise ValueError

A_EL = ELEMS[16]; B_EL = ELEMS[12]; P_EL = ELEMS[1]
assert A_EL == (((0, 1), (2, 2)), 1) and B_EL == (((0, 1), (2, 1)), 2) and P_EL == (I2, 1)
assert all(gmul(P_EL, g) == gmul(g, P_EL) for g in ELEMS), 'p must be central'

# --- the chart (row-major matrices, directions) -------------------------------------------------------------------
CHART = {
    'a': ([[0, 0, 2, 0], [0, 0, 2, 1], [1, 0, 2, 1], [2, 1, 2, 0]], [1, 1, 0, 0]),
    'A': ([[0, 2, 1, 0], [2, 2, 1, 1], [2, 0, 0, 0], [2, 1, 0, 0]], [2, 1, 2, 0]),   # a^-1
    'b': ([[2, 0, 0, 2], [1, 2, 0, 0], [2, 2, 0, 2], [2, 1, 1, 2]], [2, 0, 0, 2]),
    'B': ([[1, 1, 2, 0], [1, 0, 2, 0], [1, 0, 1, 1], [1, 2, 1, 0]], [2, 2, 1, 2]),   # b^-1
}
C_P = [2, 0, 2, 2]   # plain constant, M_p = M_{p^-1} = -I

def digits(x): return [(x // 3 ** i) % 3 for i in range(4)]
def undigits(v): return sum(v[i] * 3 ** i for i in range(4))
def matvec(A, v): return [sum(A[i][j] * v[j] for j in range(4)) % 3 for i in range(4)]
def mat_mul4(X, Y): return [[sum(X[i][k] * Y[k][j] for k in range(4)) % 3 for j in range(4)] for i in range(4)]

def check_chart():
    I4 = [[int(i == j) for j in range(4)] for i in range(4)]
    for t, ti in (('a', 'A'), ('b', 'B')):
        At, dt = CHART[t]; Ai, di = CHART[ti]
        assert mat_mul4(At, Ai) == I4, 'inverse pair'
        assert matvec(Ai, dt) == di, 'inverse direction'

def build():
    check_chart()
    SYM = [('a', A_EL), ('A', ginv(A_EL)), ('b', B_EL), ('B', ginv(B_EL))]
    edges = set()
    for g in ELEMS:
        gi = INDEX[g]
        for x in range(81):
            u = 81 * gi + x; xv = digits(x)
            for name, el in SYM:
                w = INDEX[gmul(g, el)]; A, d = CHART[name]; y = matvec(A, xv)
                for lam in range(3):
                    z = [(y[i] + lam * d[i]) % 3 for i in range(4)]
                    v = 81 * w + undigits(z)
                    edges.add((min(u, v), max(u, v)))
            for el in (P_EL, ginv(P_EL)):
                w = INDEX[gmul(g, el)]
                z = [(-xv[i] + C_P[i]) % 3 for i in range(4)]
                v = 81 * w + undigits(z)
                edges.add((min(u, v), max(u, v)))
    return sorted(edges)

if __name__ == '__main__':
    edges = build()
    n = 120 * 81
    deg = [0] * n
    for u, v in edges: deg[u] += 1; deg[v] += 1
    assert len(edges) == 7 * n and min(deg) == max(deg) == 14, (len(edges), min(deg), max(deg))
    with open(sys.argv[1], 'w') as f:
        for u, v in edges: f.write('%d %d\n' % (u, v))
    print('vertices', n, 'edges', len(edges))
