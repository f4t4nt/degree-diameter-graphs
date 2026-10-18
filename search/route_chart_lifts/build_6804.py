#!/usr/bin/env python3
"""Self-contained construction of the (13,4) graph on 6804 vertices: a route-chart lift with the five-symbol MIXED
alphabet over F_3^4 (four F_3-line symbols a, a^-1, b, b^-1 with 3 neighbours each, one plain self-inverse symbol s
with 1 neighbour) of the Cayley graph of G = Z_7 x|_2 Z_12 (order 84).  No dependency on the search code.
Writes the sorted edge list (vertices 0..6803, one "u v" per line, u < v) to the file given as argv[1].

Group G = Z_7 x|_2 Z_12: pairs (x, y), x in Z_7, y in Z_12, (x, y) (x', y') = (x + 2^y x' mod 7, y + y' mod 12)
(2 has order 3 mod 7, so the action of Z_12 on Z_7 factors through Z_3).  Element (x, y) has the index 7 y + x
(identity 0).  Centre: {(0, 3k)} = Z_4 (2^3 = 8 = 1 mod 7); the unique central involution is z = (0, 6) = 42.
Generators: a = 8 = (1, 1) of order 12, b = 28 = (0, 4) of order 3, plain symbol z = (0, 6) (central involution).
Controller C = Cay(G, {a, a^-1, b, b^-1, z}): g -> g a, g a^-1, g b, g b^-1, g z (right multiplication), a 5-regular
graph on 84 vertices; the four line symbols form a 4-regular graph on their own; z is central so the plain edges
commute with the line edges.
Fibre F_3^4, x = (x0, x1, x2, x3) encoded as x0 + 3 x1 + 9 x2 + 27 x3; lift vertex (g, x) has index 81 g + x.
Chart: for each line symbol t a matrix A_t in GL_4(F_3) (rows below, y_i = sum_j A[i][j] x_j) and a direction d_t
with A_{a^-1} = A_a^-1, d_{a^-1} = A_a^-1 d_a (same for b); for the plain symbol M_z = -I, c_z = (0, 0, 1, 0)
(M_z^2 = I and (M_z + I) c_z = 0, so the plain edges are symmetric).
Edges: (g, x) ~ (g t, A_t x + lambda d_t) for t in {a, a^-1, b, b^-1}, lambda in F_3 (12 neighbours), and
(g, x) ~ (g z, -x + c_z) (1 neighbour): 13-regular, 6804 vertices, diameter 4.
"""
import sys

# --- the group ---------------------------------------------------------------------------------------------------
M7, M12 = 7, 12
RP = [pow(2, y, M7) for y in range(M12)]
ELEMS = [(x, y) for y in range(M12) for x in range(M7)]
INDEX = {g: k for k, g in enumerate(ELEMS)}
def gmul(u, v): return ((u[0] + RP[u[1]] * v[0]) % M7, (u[1] + v[1]) % M12)
def ginv(u):
    for v in ELEMS:
        if gmul(u, v) == (0, 0): return v
    raise ValueError

A_EL = ELEMS[8]; B_EL = ELEMS[28]; Z_EL = ELEMS[42]
assert A_EL == (1, 1) and B_EL == (0, 4) and Z_EL == (0, 6)
assert all(gmul(Z_EL, g) == gmul(g, Z_EL) for g in ELEMS) and gmul(Z_EL, Z_EL) == (0, 0), 'z must be a central involution'

# --- the chart (row-major matrices, directions) -------------------------------------------------------------------
CHART = {
    'a': ([[0, 1, 0, 1], [2, 0, 2, 0], [2, 1, 1, 1], [1, 0, 1, 2]], [0, 0, 1, 1]),
    'A': ([[2, 1, 1, 0], [1, 1, 0, 1], [1, 1, 2, 0], [0, 2, 0, 2]], [1, 1, 2, 2]),   # a^-1
    'b': ([[0, 0, 2, 0], [2, 1, 2, 2], [2, 1, 2, 0], [0, 1, 2, 2]], [1, 2, 2, 2]),
    'B': ([[0, 2, 0, 1], [2, 2, 1, 1], [2, 0, 0, 0], [0, 2, 1, 0]], [0, 1, 2, 0]),   # b^-1
}
C_Z = [0, 0, 1, 0]   # plain constant, M_z = -I

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
            w = INDEX[gmul(g, Z_EL)]
            z = [(-xv[i] + C_Z[i]) % 3 for i in range(4)]
            v = 81 * w + undigits(z)
            edges.add((min(u, v), max(u, v)))
    return sorted(edges)

if __name__ == '__main__':
    edges = build()
    n = 84 * 81
    deg = [0] * n
    for u, v in edges: deg[u] += 1; deg[v] += 1
    assert len(edges) == 13 * n // 2 and min(deg) == max(deg) == 13, (len(edges), min(deg), max(deg))
    with open(sys.argv[1], 'w') as f:
        for u, v in edges: f.write('%d %d\n' % (u, v))
    print('vertices', n, 'edges', len(edges))
