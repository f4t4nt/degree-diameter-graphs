"""Finite groups as multiplication tables, built explicitly in Python.

Constructors: cyclic, abelian (direct product of cyclics), direct product,
semidirect product A x| H for an abelian normal subgroup A and any group H
(all homomorphisms H -> Aut(A), found by assigning images of generators of H),
dicyclic, permutation groups (A4, S4, A5), matrix groups over Z_p (SL(2,3)),
central quotients.  Groups are deduplicated by a pynauty certificate of the
coloured Cayley-table graph, which is a complete isomorphism invariant.
Aut(G) is read off the same graph (with the identity fixed).

Element 0 is always the identity.  mul[a, b] = a*b.  The search uses right
multiplication by voltages, so left multiplication by G gives the lift its
automorphisms.
"""
import itertools
import math

import numpy as np
import pynauty

GNU = {  # number of groups of order n (OEIS A000001), used only to report coverage
    1: 1, 2: 1, 3: 1, 4: 2, 5: 1, 6: 2, 7: 1, 8: 5, 9: 2, 10: 2, 11: 1, 12: 5, 13: 1, 14: 2,
    15: 1, 16: 14, 17: 1, 18: 5, 19: 1, 20: 5, 21: 2, 22: 2, 23: 1, 24: 15, 25: 2, 26: 2,
    27: 5, 28: 4, 29: 1, 30: 4, 31: 1, 32: 51, 33: 1, 34: 2, 35: 1, 36: 14, 37: 1, 38: 2,
    39: 2, 40: 14, 41: 1, 42: 6, 43: 1, 44: 4, 45: 2, 46: 2, 47: 1, 48: 52, 49: 2, 50: 5,
    51: 1, 52: 5, 53: 1, 54: 15, 55: 2, 56: 13, 57: 2, 58: 2, 59: 1, 60: 13, 61: 1, 62: 2,
    63: 4, 64: 267, 65: 1, 66: 4, 67: 1, 68: 5, 69: 1, 70: 4, 71: 1, 72: 50, 73: 1, 74: 2,
    75: 3, 76: 4, 77: 1, 78: 6, 79: 1, 80: 52, 81: 15, 82: 2, 83: 1, 84: 15, 85: 1, 86: 2,
    87: 1, 88: 12, 89: 1, 90: 10, 91: 1, 92: 4, 93: 2, 94: 2, 95: 1, 96: 231, 97: 1, 98: 5,
    99: 2, 100: 16, 101: 1, 102: 4, 103: 1, 104: 14, 105: 2, 106: 2, 107: 1, 108: 45,
    134: 2, 136: 15, 138: 4, 140: 11, 144: 197, 198: 10, 200: 52, 204: 12, 208: 51,
    210: 12, 216: 177,
}


class Group:
    def __init__(self, name, mul):
        self.name = name
        self.mul = np.asarray(mul, dtype=np.int32)
        self.n = self.mul.shape[0]
        assert self.mul.shape == (self.n, self.n)
        assert (self.mul[0] == np.arange(self.n)).all() and (self.mul[:, 0] == np.arange(self.n)).all()
        inv = np.full(self.n, -1, dtype=np.int32)
        for a in range(self.n):
            z = np.nonzero(self.mul[a] == 0)[0]
            assert len(z) == 1
            inv[a] = z[0]
        self.inv = inv
        self._cert = None
        self._aut = None

    def check_associative(self, samples=2000):
        rng = np.random.default_rng(0)
        a = rng.integers(0, self.n, samples)
        b = rng.integers(0, self.n, samples)
        c = rng.integers(0, self.n, samples)
        return (self.mul[self.mul[a, b], c] == self.mul[a, self.mul[b, c]]).all()

    def order_of(self, a):
        k, x = 1, a
        while x != 0:
            x = self.mul[x, a]
            k += 1
        return k

    def element_orders(self):
        return [self.order_of(a) for a in range(self.n)]

    def involutions(self):
        return [a for a in range(self.n) if a != 0 and self.mul[a, a] == 0]

    def is_abelian(self):
        return (self.mul == self.mul.T).all()

    def center_size(self):
        return int(sum((self.mul[a] == self.mul[:, a]).all() for a in range(self.n)))

    def table_graph(self, fix_identity):
        """Coloured graph encoding the Cayley table: vertices r_a, c_b, p_c (3n)
        and one vertex per cell (a,b) adjacent to r_a, c_b, p_{ab}."""
        n = self.n
        adj = {}
        for a in range(n):
            for b in range(n):
                cell = 3 * n + a * n + b
                adj[cell] = [a, n + b, 2 * n + int(self.mul[a, b])]
        g = pynauty.Graph(3 * n + n * n, adjacency_dict=adj)
        if fix_identity:
            colouring = [{0}, {n}, {2 * n}, set(range(1, n)), set(range(n + 1, 2 * n)),
                         set(range(2 * n + 1, 3 * n)), set(range(3 * n, 3 * n + n * n))]
        else:
            colouring = [set(range(n)), set(range(n, 2 * n)), set(range(2 * n, 3 * n)),
                         set(range(3 * n, 3 * n + n * n))]
        g.set_vertex_coloring(colouring)
        return g

    def certificate(self):
        if self._cert is None:
            self._cert = pynauty.certificate(self.table_graph(False))
        return self._cert

    def element_signature(self):
        """Per-element isomorphism-invariant signature (order, #square roots,
        conjugacy class size, order of the square)."""
        if getattr(self, "_esig", None) is not None:
            return self._esig
        n = self.n
        orders = self.element_orders()
        sq = self.mul[np.arange(n), np.arange(n)]
        nsqrt = np.bincount(sq, minlength=n)
        g = np.arange(n)[:, None]
        x = np.arange(n)[None, :]
        conj = self.mul[self.mul[g, x], self.inv[g]]  # conj[g, x] = g x g^-1
        ccsize = [len(set(conj[:, i].tolist())) for i in range(n)]
        self._esig = [(orders[i], int(nsqrt[i]), ccsize[i], orders[int(sq[i])]) for i in range(n)]
        return self._esig

    def derived_subgroup_size(self):
        n = self.n
        g = np.arange(n)[:, None]
        x = np.arange(n)[None, :]
        comm = self.mul[self.mul[g, x], self.mul[self.inv[g], self.inv[x]]]
        gens = sorted(set(comm.reshape(-1).tolist()))
        return len(closure(gens, lambda a, b: int(self.mul[a, b]), 0))

    def signature(self):
        if getattr(self, "_sig", None) is None:
            from collections import Counter
            es = Counter(self.element_signature())
            self._sig = (self.n, tuple(sorted(es.items())), self.derived_subgroup_size(), self.center_size())
        return self._sig

    def aut_generators(self):
        """Generators of Aut(G) as permutations of {0..n-1}."""
        gens, _, _, _, _ = pynauty.autgrp(self.table_graph(True))
        return [tuple(g[: self.n]) for g in gens]

    def aut_group(self, limit=60000):
        """Aut(G) expanded from generators (or a subgroup if bigger than limit)."""
        if self._aut is None:
            gens = self.aut_generators()
            self._aut = expand_perm_group(gens, self.n, limit)
        return self._aut

    def aut_order(self):
        _, a, b, _, _ = pynauty.autgrp(self.table_graph(True))
        return int(round(a * 10**b))

def expand_perm_group(gens, n, limit):
    ident = tuple(range(n))
    seen = {ident}
    frontier = [ident]
    gens = [tuple(g) for g in gens]
    while frontier:
        nxt = []
        for g in frontier:
            for s in gens:
                h = tuple(s[g[x]] for x in range(n))
                if h not in seen:
                    seen.add(h)
                    nxt.append(h)
                    if len(seen) >= limit:
                        return sorted(seen)
        frontier = nxt
    return sorted(seen)


# ---------------------------------------------------------------- constructors

def from_elements(name, elems, op):
    """Group from a list of hashable elements (elems[0] = identity) and a binary op."""
    idx = {e: i for i, e in enumerate(elems)}
    n = len(elems)
    mul = np.zeros((n, n), dtype=np.int32)
    for i, a in enumerate(elems):
        for j, b in enumerate(elems):
            mul[i, j] = idx[op(a, b)]
    return Group(name, mul)


def closure(gens, op, ident):
    elems = [ident]
    seen = {ident}
    frontier = [ident]
    while frontier:
        nxt = []
        for a in frontier:
            for g in gens:
                b = op(a, g)
                if b not in seen:
                    seen.add(b)
                    elems.append(b)
                    nxt.append(b)
        frontier = nxt
    return elems


def cyclic(n):
    return from_elements(f"Z{n}", list(range(n)), lambda a, b: (a + b) % n)


def abelian(ms):
    ms = tuple(ms)
    if len(ms) == 1:
        return cyclic(ms[0])
    elems = list(itertools.product(*[range(m) for m in ms]))
    return from_elements("x".join(f"Z{m}" for m in ms), elems,
                         lambda a, b: tuple((x + y) % m for x, y, m in zip(a, b, ms)))


def direct_product(G, H):
    n, m = G.n, H.n
    mul = np.zeros((n * m, n * m), dtype=np.int32)
    gm = G.mul
    hm = H.mul
    for a in range(n):
        for b in range(m):
            # element index a*m + b
            mul[a * m + b] = (gm[a][:, None] * m + hm[b][None, :]).reshape(-1)
    return Group(f"({G.name})x({H.name})", mul)


def perm_group(name, gens):
    """Group generated by permutations (tuples)."""
    deg = len(gens[0])
    ident = tuple(range(deg))
    op = lambda p, q: tuple(p[q[i]] for i in range(deg))
    elems = closure(gens, op, ident)
    return from_elements(name, elems, op)


def matrix_group(name, gens, p):
    gens = [tuple(tuple(int(x) % p for x in row) for row in g) for g in gens]
    d = len(gens[0])
    ident = tuple(tuple(1 if i == j else 0 for j in range(d)) for i in range(d))

    def op(a, b):
        return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(d)) % p for j in range(d)) for i in range(d))

    elems = closure(gens, op, ident)
    return from_elements(name, elems, op)


def dicyclic(m):
    """Dic_m of order 4m: <a, x | a^{2m}=1, x^2=a^m, x a x^-1 = a^-1>; elements a^i x^j."""
    n2 = 2 * m
    elems = [(i, j) for j in range(2) for i in range(n2)]

    def op(u, v):
        i, j = u
        k, l = v
        kk = k if j == 0 else (-k) % n2
        s = (i + kk) % n2
        if j == 1 and l == 1:
            return ((s + m) % n2, 0)
        return (s, (j + l) % 2)

    return from_elements(f"Dic{m}", elems, op)


def semidirect(A, H, phi, name=None):
    """A x| H where phi[h] is a permutation of A (an automorphism) for each h in H,
    with phi a homomorphism H -> Aut(A).  (a,h)(b,k) = (a * phi_h(b), h k)."""
    n, m = A.n, H.n
    phi = np.asarray(phi, dtype=np.int32)
    mul = np.zeros((n * m, n * m), dtype=np.int32)
    for a in range(n):
        for h in range(m):
            # (a,h)(b,k) for all b,k
            pb = phi[h]  # phi_h(b) for all b
            mul[a * m + h] = (A.mul[a][pb][:, None] * m + H.mul[h][None, :]).reshape(-1)
    return Group(name or f"({A.name})x|({H.name})", mul)


def generators_of(G):
    """A small generating set (greedy)."""
    gens = []
    sub = {0}
    orders = G.element_orders()
    cand = sorted(range(1, G.n), key=lambda a: -orders[a])
    for a in cand:
        if a in sub:
            continue
        gens.append(a)
        sub = set(closure(gens, lambda x, y: int(G.mul[x, y]), 0))
        if len(sub) == G.n:
            break
    return gens


def homomorphisms(H, target_perms, n_target_elems):
    """All homomorphisms H -> T where T is a group of permutations (list of tuples),
    returned as arrays phi[h] = permutation.  Enumerated by images of a generating
    set of H, verified by extending along the Cayley graph of H."""
    gens = generators_of(H)
    perm_index = {p: i for i, p in enumerate(target_perms)}
    ident = tuple(range(n_target_elems))
    # words: express each element of H as (prev element, generator) via BFS
    parent = {0: None}
    order = [0]
    frontier = [0]
    while frontier:
        nxt = []
        for a in frontier:
            for gi, g in enumerate(gens):
                b = int(H.mul[a, g])
                if b not in parent:
                    parent[b] = (a, gi)
                    order.append(b)
                    nxt.append(b)
        frontier = nxt
    assert len(order) == H.n
    homs = []
    for images in itertools.product(range(len(target_perms)), repeat=len(gens)):
        img = [None] * H.n
        img[0] = ident
        for b in order[1:]:
            a, gi = parent[b]
            pa = img[a]
            pg = target_perms[images[gi]]
            img[b] = tuple(pa[pg[x]] for x in range(n_target_elems))  # phi(a g) = phi(a) o phi(g)
        # verify homomorphism on all pairs with generators (sufficient)
        ok = True
        for a in range(H.n):
            for gi, g in enumerate(gens):
                b = int(H.mul[a, g])
                pa = img[a]
                pg = target_perms[images[gi]]
                if img[b] != tuple(pa[pg[x]] for x in range(n_target_elems)):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            homs.append(np.array(img, dtype=np.int32))
    return homs


def isomorphic(G, H):
    """Exact isomorphism test: map a generating set of G to elements of H with the
    same element signature and check the extension is a bijective homomorphism."""
    if G.n != H.n or G.signature() != H.signature():
        return False
    gens = generators_of(G)
    sigG = G.element_signature()
    sigH = H.element_signature()
    cands = [[h for h in range(H.n) if sigH[h] == sigG[g]] for g in gens]
    # BFS words for G
    parent = {0: None}
    order = [0]
    frontier = [0]
    while frontier:
        nxt = []
        for a in frontier:
            for gi, g in enumerate(gens):
                b = int(G.mul[a, g])
                if b not in parent:
                    parent[b] = (a, gi)
                    order.append(b)
                    nxt.append(b)
        frontier = nxt
    Gm, Hm = G.mul, H.mul
    for images in itertools.product(*cands):
        img = [-1] * G.n
        img[0] = 0
        for b in order[1:]:
            a, gi = parent[b]
            img[b] = int(Hm[img[a], images[gi]])
        if len(set(img)) != G.n:
            continue
        ok = True
        for a in range(G.n):
            for gi, g in enumerate(gens):
                if img[int(Gm[a, g])] != Hm[img[a], images[gi]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return True
    return False


def quotient(G, N_elems, name):
    """G / N for a normal subgroup N given by its elements."""
    N = sorted(N_elems)
    coset_of = {}
    reps = []
    for a in range(G.n):
        if a in coset_of:
            continue
        c = len(reps)
        reps.append(a)
        for x in N:
            coset_of[int(G.mul[a, x])] = c
    k = len(reps)
    mul = np.zeros((k, k), dtype=np.int32)
    for i, a in enumerate(reps):
        for j, b in enumerate(reps):
            mul[i, j] = coset_of[int(G.mul[a, b])]
    # make identity coset index 0 (it is, since a=0 is processed first)
    return Group(name, mul)


# ------------------------------------------------------------- group pools

def abelian_types(n):
    """All abelian groups of order n as tuples of cyclic factors (primary decomposition)."""
    fac = factorize(n)
    per_prime = []
    for p, e in fac.items():
        per_prime.append([[p**part for part in parts] for parts in partitions(e)])
    out = []
    for combo in itertools.product(*per_prime):
        out.append(tuple(sorted(itertools.chain.from_iterable(combo))))
    return out


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield []
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield [k] + rest


def factorize(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


SPECIAL = {}


def special_groups():
    if SPECIAL:
        return SPECIAL
    SPECIAL[12] = [perm_group("A4", [(1, 2, 0, 3), (0, 2, 3, 1)])]
    SPECIAL[24] = [perm_group("S4", [(1, 0, 2, 3), (1, 2, 3, 0)]),
                   matrix_group("SL(2,3)", [((1, 1), (0, 1)), ((0, -1), (1, 0))], 3)]
    SPECIAL[60] = [perm_group("A5", [(1, 2, 0, 3, 4), (0, 1, 3, 4, 2), (1, 0, 3, 2, 4)])]
    SPECIAL[8] = [dicyclic(2)]
    SPECIAL[48] = [matrix_group("GL(2,3)", [((1, 1), (0, 1)), ((0, -1), (1, 0)), ((2, 0), (0, 1))], 3)]
    SPECIAL[20] = [perm_group("F20", [(1, 2, 3, 4, 0), (0, 2, 4, 1, 3)])]
    SPECIAL[21] = [perm_group("F21", [(1, 2, 3, 4, 5, 6, 0), (0, 2, 4, 6, 1, 3, 5)])]
    SPECIAL[18] = [perm_group("Z3^2x|Z2", [(1, 2, 0, 3, 4, 5), (0, 1, 2, 4, 5, 3), (0, 2, 1, 3, 5, 4)])]
    return SPECIAL


_POOL = {}


class Dedup:
    """Isomorphism-class collector: cheap signature first, exact test on collision."""

    def __init__(self):
        self.by_sig = {}
        self.groups = []

    def add(self, G):
        s = G.signature()
        bucket = self.by_sig.setdefault(s, [])
        for H in bucket:
            if isomorphic(H, G):
                return False
        bucket.append(G)
        self.groups.append(G)
        return True


def pool(n, max_aut=60000, verbose=False, cyclic_complements_only=None, time_budget=None):
    """Non-isomorphic groups of order n constructible here."""
    if n in _POOL:
        return _POOL[n]
    import time
    t0 = time.time()
    if cyclic_complements_only is None:
        cyclic_complements_only = n > 110
    dd = Dedup()
    for t in abelian_types(n):
        dd.add(abelian(t))
    if n % 4 == 0 and n >= 8:
        dd.add(dicyclic(n // 4))
    for g in special_groups().get(n, []):
        dd.add(g)
    # direct products of smaller pool members
    for d in divisors(n):
        if 1 < d < n and d <= n // d:
            for G in pool(d, max_aut, verbose):
                for H in pool(n // d, max_aut, verbose):
                    dd.add(direct_product(G, H))
    # semidirect products A x| H, A abelian
    for d in divisors(n):
        if 1 < d < n:
            for t in abelian_types(d):
                A = abelian(t)
                try:
                    aut = A.aut_group(max_aut)
                except MemoryError:
                    continue
                if len(aut) >= max_aut:
                    continue
                for H in pool(n // d, max_aut, verbose):
                    gens = generators_of(H)
                    if cyclic_complements_only and len(gens) > 1:
                        continue
                    if len(aut) ** len(gens) > 300_000:
                        continue
                    for phi in homomorphisms(H, aut, A.n):
                        dd.add(semidirect(A, H, phi))
                    if time_budget and time.time() - t0 > time_budget:
                        break
    out = dd.groups
    out.sort(key=lambda G: (not G.is_abelian(), G.name))
    for i, G in enumerate(out):
        G.index = i
    _POOL[n] = out
    if verbose:
        print(f"order {n}: {len(out)} groups found (known total {GNU.get(n, '?')}) in {time.time() - t0:.1f}s", flush=True)
    return out


def describe(G, with_aut=True):
    from collections import Counter
    orders = Counter(G.element_orders())
    a = f"|Aut|={G.aut_order()}, " if with_aut else ""
    return (f"{G.name}: order {G.n}, abelian={G.is_abelian()}, |Z|={G.center_size()}, "
            f"{a}involutions={len(G.involutions())}, "
            f"element orders={dict(sorted(orders.items()))}")


def save_group(G, path, with_aut=True):
    """Group file for lift.c: n, table, then the automorphisms (or just the identity)."""
    aut = [tuple(range(G.n))]
    if with_aut:
        try:
            aut = G.aut_group()
        except MemoryError:
            pass
    with open(path, "w") as f:
        f.write(f"{G.n}\n")
        for row in G.mul:
            f.write(" ".join(map(str, row)) + "\n")
        f.write(f"{len(aut)}\n")
        for p in aut:
            f.write(" ".join(map(str, p)) + "\n")
    return len(aut)


if __name__ == "__main__":
    import os
    import sys
    outdir = sys.argv[1]
    with_aut = os.environ.get("NO_AUT") is None
    for n in map(int, sys.argv[2:]):
        gs = pool(n, verbose=True)
        d = os.path.join(outdir, f"order_{n}")
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "index.txt"), "w") as idx:
            idx.write(f"# order {n}: {len(gs)} groups found; known number of groups {GNU.get(n, '?')}\n")
            for i, G in enumerate(gs):
                naut = save_group(G, os.path.join(d, f"g{i}.txt"), with_aut and n <= 110)
                line = f"g{i} naut_used={naut} {describe(G, with_aut=False)}"
                idx.write(line + "\n")
                print("  ", line, flush=True)
