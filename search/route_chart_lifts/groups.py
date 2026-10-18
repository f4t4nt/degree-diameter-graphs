#!/usr/bin/env python3
"""Small groups as multiplication tables (text: n, then n*n entries T[g*n+h] = g*h, identity 0).
catalogue(order) lists constructible groups of that order: cyclic, abelian, dihedral, dicyclic, metacyclic
Z_m x|_r Z_M (all r up to automorphism), direct products with cyclic groups, and named exceptions (A4, S4, SL(2,3),
Heis(p), Z_p^2 x| Z_p, Z_3^2 x| Z_4, A4 x Z_n, S4 x Z_2, ...).  Complete for all orders <= 31 except possibly a few
2-groups; for larger orders it is the metacyclic + abelian + product subfamily.
usage: groups.py write ORDER OUTDIR   (writes OUTDIR/<name>.tbl, prints names)"""
import sys, os, itertools
from math import gcd

class Group:
    def __init__(self, name, elems, mul, e):
        self.name = name; self.elems = list(elems); self.mul = mul; self.e = e
        if self.elems[0] != e:
            self.elems.remove(e); self.elems.insert(0, e)
        self.index = {x: i for i, x in enumerate(self.elems)}
        self.n = len(self.elems)
    def table(self):
        idx = self.index; el = self.elems
        return [[idx[self.mul(g, h)] for h in el] for g in el]
    def write(self, fn):
        T = self.table()
        with open(fn, 'w') as f:
            f.write('%d\n' % self.n)
            for row in T: f.write(' '.join(map(str, row)) + '\n')
    def invariants(self):
        T = self.table(); n = self.n
        inv = [row.index(0) for row in T]
        ninv = sum(1 for g in range(1, n) if inv[g] == g)
        centre = sum(1 for g in range(n) if all(T[g][h] == T[h][g] for h in range(n)))
        abelian = centre == n
        # exponent / element order histogram
        orders = {}
        for g in range(n):
            x, o = g, 1
            while x != 0: x = T[x][g]; o += 1
            orders[o] = orders.get(o, 0) + 1
        cen_orders = {}
        for g in range(n):
            if all(T[g][h] == T[h][g] for h in range(n)):
                x, o = g, 1
                while x != 0: x = T[x][g]; o += 1
                cen_orders[o] = cen_orders.get(o, 0) + 1
        return ninv, centre, abelian, tuple(sorted(orders.items())), tuple(sorted(cen_orders.items()))

def perm_group(name, gens):
    deg = len(gens[0]); e = tuple(range(deg))
    elems = {e}; frontier = [e]
    while frontier:
        new = []
        for p in frontier:
            for g in gens:
                q = tuple(p[g[i]] for i in range(deg))
                if q not in elems: elems.add(q); new.append(q)
        frontier = new
    return Group(name, sorted(elems), lambda p, q: tuple(p[q[i]] for i in range(deg)), e)

def cyclic(n): return Group('Z%d' % n, range(n), lambda a, b: (a + b) % n, 0)

def direct(G, H, name=None):
    return Group(name or '%sx%s' % (G.name, H.name), [(g, h) for g in G.elems for h in H.elems],
                 lambda a, b: (G.mul(a[0], b[0]), H.mul(a[1], b[1])), (G.e, H.e))

def abelian(orders):
    G = cyclic(orders[0])
    for o in orders[1:]: G = direct(G, cyclic(o))
    G.name = 'Z' + 'xZ'.join(map(str, orders))
    return G

def metacyclic(m, M, r, name=None):
    assert pow(r, M, m) == 1
    rp = [pow(r, y, m) for y in range(M)]
    return Group(name or 'Z%d:%d:Z%d' % (m, r, M), [(x, y) for y in range(M) for x in range(m)],
                 lambda a, b: ((a[0] + rp[a[1]] * b[0]) % m, (a[1] + b[1]) % M), (0, 0))

def dihedral(n): return metacyclic(n, 2, n - 1, 'D%d' % n)

def dicyclic(n):
    """Dic_n of order 4n: <a, x | a^2n = 1, x^2 = a^n, x a x^-1 = a^-1>; elements a^i x^j, j in {0,1}."""
    m = 2 * n
    def mul(u, v):
        i, j = u; k, l = v
        if j == 0: return ((i + k) % m, l)
        # x a^k = a^-k x ; x x = a^n
        if l == 0: return ((i - k) % m, 1)
        return ((i - k + n) % m, 0)
    return Group('Dic%d' % n, [(i, j) for j in range(2) for i in range(m)], mul, (0, 0))

def gen_dihedral(A, name):
    """A x| Z_2 with inversion (A abelian group given as Group with elems tuples/ints)."""
    def mul(u, v):
        a, s = u; b, t = v
        if s == 0: return (A.mul(a, b), t)
        return (A.mul(a, A_inv(A, b)), (s + t) % 2)
    return Group(name, [(a, s) for s in range(2) for a in A.elems], mul, (A.e, 0))

def A_inv(A, b):
    for x in A.elems:
        if A.mul(b, x) == A.e: return x

def heis(p):
    return Group('Heis%d' % p, [(a, b, c) for a in range(p) for b in range(p) for c in range(p)],
                 lambda u, v: ((u[0] + v[0]) % p, (u[1] + v[1]) % p, (u[2] + v[2] + u[0] * v[1]) % p), (0, 0, 0))

def zp2_zp(p):
    """Z_{p^2} x| Z_p, generator acting by 1 + p."""
    return metacyclic(p * p, p, 1 + p, 'Z%d:Z%d' % (p * p, p))

def A4(): return perm_group('A4', [(1, 2, 0, 3), (1, 0, 3, 2)])
def S4(): return perm_group('S4', [(1, 0, 2, 3), (1, 2, 3, 0)])
def Q8(): return dicyclic(2)
def SL23():
    # SL(2,3) as 2x2 matrices over F_3
    els = [((a, b), (c, d)) for a in range(3) for b in range(3) for c in range(3) for d in range(3) if (a * d - b * c) % 3 == 1]
    def mul(X, Y):
        (a, b), (c, d) = X; (e, f), (g, h) = Y
        return (((a * e + b * g) % 3, (a * f + b * h) % 3), ((c * e + d * g) % 3, (c * f + d * h) % 3))
    return Group('SL23', els, mul, ((1, 0), (0, 1)))

def z3sq_z4():
    """Z_3^2 x| Z_4 (order 36), generator acting by the matrix [[0,-1],[1,0]] of order 4."""
    def act(v, j):
        for _ in range(j): v = ((-v[1]) % 3, v[0])
        return v
    return Group('Z3^2:Z4', [((x, y), j) for j in range(4) for x in range(3) for y in range(3)],
                 lambda u, v: (tuple((a + b) % 3 for a, b in zip(u[0], act(v[0], u[1]))), (u[1] + v[1]) % 4), ((0, 0), 0))

def z3sq_z2():  # generalized dihedral of Z_3^2 (order 18)
    return gen_dihedral(abelian([3, 3]), 'Z3^2:Z2')

def unit_orders(m):
    """representatives r of units mod m up to the equivalence r ~ r^j (gcd(j, ord r) = 1), with ord r > 1."""
    reps = {}
    for r in range(2, m):
        if gcd(r, m) != 1: continue
        o, x = 1, r
        while x != 1: x = (x * r) % m; o += 1
        # canonical: the cyclic subgroup generated by r
        sub = frozenset(pow(r, j, m) for j in range(o))
        if sub not in reps: reps[sub] = (r, o)
    return list(reps.values())

def dedupe(groups):
    seen = {}; out = []
    for G in groups:
        if G.n == 0: continue
        key = G.invariants()
        if key in seen: continue
        seen[key] = G.name; out.append(G)
    return out

def catalogue(n):
    gs = [cyclic(n)]
    # abelian: all partitions of the prime factorisation into cyclic factors (via divisor chains)
    def abelian_types(n):
        res = set()
        def rec(rem, parts):
            if rem == 1: res.add(tuple(sorted(parts, reverse=True))); return
            for d in range(2, rem + 1):
                if rem % d == 0 and (not parts or True): rec(rem // d, parts + [d])
        rec(n, [])
        return res
    for t in abelian_types(n):
        if len(t) > 1: gs.append(abelian(list(t)))
    if n % 2 == 0 and n >= 6: gs.append(dihedral(n // 2))
    if n % 4 == 0 and n >= 8: gs.append(dicyclic(n // 4))
    # metacyclic Z_m x|_r Z_M for all factorisations n = m M and all r classes
    for m in range(3, n):
        if n % m: continue
        M = n // m
        if M < 2: continue
        for r, o in unit_orders(m):
            if M % o == 0: gs.append(metacyclic(m, M, r))
    # products of nonabelian small groups with cyclic groups
    specials = {12: [A4()], 24: [S4(), SL23(), direct(A4(), cyclic(2))], 8: [Q8()], 27: [heis(3), zp2_zp(3)],
                18: [z3sq_z2()], 36: [z3sq_z4(), direct(A4(), cyclic(3))],
                48: [direct(S4(), cyclic(2)), direct(SL23(), cyclic(2)), direct(A4(), cyclic(4)), direct(A4(), abelian([2, 2]))],
                16: [direct(Q8(), cyclic(2)), direct(dihedral(4), cyclic(2))], 32: [direct(Q8(), cyclic(4)), direct(dihedral(4), cyclic(4)), direct(Q8(), abelian([2, 2])), direct(dihedral(4), abelian([2, 2]))],
                54: [direct(heis(3), cyclic(2)), direct(zp2_zp(3), cyclic(2)), gen_dihedral(abelian([3, 3, 3]), 'Z3^3:Z2'), gen_dihedral(abelian([9, 3]), 'Z9xZ3:Z2')],
                60: [perm_group('A5', [(1, 2, 0, 3, 4), (0, 1, 3, 4, 2)]), direct(A4(), cyclic(5))],
                72: [direct(A4(), cyclic(6)), direct(S4(), cyclic(3)), direct(SL23(), cyclic(3)), direct(z3sq_z4(), cyclic(2)), direct(z3sq_z2(), cyclic(4))],
                80: [direct(A4(), abelian([5])), ],
                81: [heis(3) and direct(heis(3), cyclic(3)), direct(zp2_zp(3), cyclic(3))],
                96: [direct(S4(), cyclic(4)), direct(SL23(), cyclic(4)), direct(A4(), cyclic(8)), direct(S4(), abelian([2, 2]))],
                100: [direct(dicyclic(5), cyclic(5)), direct(dihedral(5), cyclic(10)), direct(metacyclic(5, 4, 2), cyclic(5))],
                108: [direct(heis(3), cyclic(4)), direct(zp2_zp(3), cyclic(4)), direct(A4(), cyclic(9)), direct(A4(), abelian([3, 3])), direct(z3sq_z4(), cyclic(3))]}
    if n in specials: gs += specials[n]
    # products of dihedral / dicyclic / metacyclic of smaller order with cyclic groups
    for d in range(2, n):
        if n % d: continue
        c = n // d
        if c == 1: continue
        for H in catalogue_nonab_small(d):
            for t in abelian_types(c):
                gs.append(direct(H, abelian(list(t))))
    return dedupe(gs)

_small_cache = {}
def catalogue_nonab_small(d):
    if d in _small_cache: return _small_cache[d]
    gs = []
    if d % 2 == 0 and d >= 6: gs.append(dihedral(d // 2))
    if d % 4 == 0 and d >= 8: gs.append(dicyclic(d // 4))
    for m in range(3, d):
        if d % m: continue
        M = d // m
        if M < 2: continue
        for r, o in unit_orders(m):
            if M % o == 0: gs.append(metacyclic(m, M, r))
    if d == 12: gs.append(A4())
    if d == 24: gs += [S4(), SL23(), direct(A4(), cyclic(2))]
    if d == 27: gs += [heis(3), zp2_zp(3)]
    if d == 18: gs.append(z3sq_z2())
    if d == 36: gs += [z3sq_z4(), direct(A4(), cyclic(3))]
    if d == 16: gs += [direct(Q8(), cyclic(2)), direct(dihedral(4), cyclic(2))]
    if d == 48: gs += [direct(S4(), cyclic(2)), direct(SL23(), cyclic(2)), direct(A4(), cyclic(4)), direct(A4(), abelian([2, 2]))]
    if d == 54: gs += [direct(heis(3), cyclic(2)), direct(zp2_zp(3), cyclic(2)), gen_dihedral(abelian([3, 3, 3]), 'Z3^3:Z2'), gen_dihedral(abelian([9, 3]), 'Z9xZ3:Z2')]
    if d == 72: gs += [direct(SL23(), cyclic(3)), direct(S4(), cyclic(3)), direct(z3sq_z4(), cyclic(2))]
    if d == 60: gs.append(perm_group('A5', [(1, 2, 0, 3, 4), (0, 1, 3, 4, 2)]))
    _small_cache[d] = dedupe(gs)
    return _small_cache[d]

if __name__ == '__main__':
    if sys.argv[1] == 'write':
        n = int(sys.argv[2]); outdir = sys.argv[3]
        os.makedirs(outdir, exist_ok=True)
        for G in catalogue(n):
            fn = os.path.join(outdir, '%d_%s.tbl' % (n, G.name.replace('/', '_')))
            G.write(fn)
            ninv, centre, ab, orders, _ = G.invariants()
            print('%s inv=%d centre=%d abelian=%s' % (fn, ninv, centre, ab))
