#!/usr/bin/env python3
"""enum_family.py: the groups G = K x| C_M, |G| in (LO, HI], K odd = a direct product of blocks, alpha = a product of block
automorphisms (one conjugacy-class representative per block automorphism class), faithful (ord alpha = M) or with kernel
exactly <y^(M/2)> (one central involution), cheap filters (involutions <= IMAX, no alpha-fixed centre of K beyond the
allowed kernel), then the full numpy build (gf.py) for |Z(G)| and the final list.

Blocks for a prime power p^e of |K| (the p-part is one block; frob blocks span p^2 x k):
  e=1: Z_p;  e=2: Z_p^2, F_p^2;  e=3: Z_p^3, Z_p^2 x Z_p, F_p^3 (p<=7), p^{1+2}_+ (p<=7);  e=4: Z_p^4, Z_p^3 x Z_p, Z_p^2 x Z_p^2,
  F_p^4 (p=3, squarefree characteristic polynomial only);  frob: F_p^2 x| <B>, B fixed-point-free of odd order k >= 3,
  k | p^2-1 (one B per k), automorphisms (A, e) with A in N(<B>), plus translation twists when A has eigenvalue 1.
Block automorphism classes: cyc: units; F_p^e: GL(e,p) conjugacy classes (rational canonical forms); frob: (A up to
N(<B>)-conjugacy, e) + twists; heis: GL(2,p) classes (phi_A) + twists; ab2: brute force Aut and conjugacy.
Reduction: alpha ~ alpha^j (j unit mod M) when every block key is exact (untwisted); twisted variants are kept.

usage: enum_family.py LO HI [--mmin 4] [--mmax 150] [--imax 31] [--zmax 2] [--out specs.jsonl] [--count-only]
"""
import itertools
import json
import sys
from math import gcd
from collections import defaultdict

import numpy as np

import gf as grp  # noqa: E402  (gf.py: the group-file builder; 'grp' is a stdlib module name)


def factor(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def lcm(a, b):
    return a * b // gcd(a, b)


def mult_order(u, q):
    if q == 1:
        return 1
    o, x = 1, u % q
    while x != 1:
        x = x * u % q; o += 1
    return o


# ---------------- polynomials over F_p (coefficient lists, low degree first) ----------------
def poly_mul(a, b, p):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] = (r[i + j] + x * y) % p
    return r


def charpoly(A, p):
    """det(xI - A) over F_p by the permutation expansion (n <= 4: at most 24 terms), low-degree-first coefficients"""
    n = len(A)
    import itertools as it
    res = [0] * (n + 1)
    for perm in it.permutations(range(n)):
        # sign
        sgn = 1
        pr = list(perm)
        for a_ in range(n):
            for b_ in range(a_ + 1, n):
                if pr[a_] > pr[b_]:
                    sgn = -sgn
        term = [1]
        for r in range(n):
            c = perm[r]
            # entry of xI - A at (r, c): x*[r==c] - A[r][c]
            term = poly_mul(term, [(-A[r][c]) % p, 1 if r == c else 0], p)
        for d_, v in enumerate(term):
            res[d_] = (res[d_] + sgn * v) % p
    return res


def poly_eval_matrix(poly, A, p):
    n = len(A)
    R = [[0] * n for _ in range(n)]
    P = [[int(i == j) for j in range(n)] for i in range(n)]
    for c in poly:
        for i in range(n):
            for j in range(n):
                R[i][j] = (R[i][j] + c * P[i][j]) % p
        P = grp.matmul(P, A, p)
    return R


def is_zero(Mx):
    return all(v == 0 for row in Mx for v in row)


def rank_mod_p(Mx, p):
    Mx = [row[:] for row in Mx]
    n, m = len(Mx), len(Mx[0])
    r = 0
    for c in range(m):
        piv = next((i for i in range(r, n) if Mx[i][c] % p), None)
        if piv is None:
            continue
        Mx[r], Mx[piv] = Mx[piv], Mx[r]
        invp = pow(Mx[r][c], -1, p)
        Mx[r] = [v * invp % p for v in Mx[r]]
        for i in range(n):
            if i != r and Mx[i][c] % p:
                f = Mx[i][c]
                Mx[i] = [(a - f * b) % p for a, b in zip(Mx[i], Mx[r])]
        r += 1
    return r


def monic_polys(deg, p, nonzero_const=True):
    for coeffs in itertools.product(range(p), repeat=deg):
        if nonzero_const and coeffs[0] == 0:
            continue
        yield list(coeffs) + [1]


def poly_divides(a, b, p):
    """a | b over F_p (both low-first)"""
    b = b[:]
    while len(b) >= len(a) and any(b):
        if b[-1] == 0:
            b.pop(); continue
        f = b[-1] * pow(a[-1], -1, p) % p
        sh = len(b) - len(a)
        for i, c in enumerate(a):
            b[sh + i] = (b[sh + i] - f * c) % p
        b.pop()
    return not any(b)


def minpoly_key(A, p):
    """(charpoly, minpoly) with minpoly found as the smallest-degree monic divisor of charpoly annihilating A"""
    cp = charpoly(A, p)
    n = len(A)
    for d in range(1, n + 1):
        for f in monic_polys(d, p, nonzero_const=False):
            if poly_divides(f, cp, p) and is_zero(poly_eval_matrix(f, A, p)):
                return tuple(cp), tuple(f)
    raise RuntimeError


def gl_key(A, p):
    """conjugacy class key in GL(e, p), e <= 3 exact ((charpoly, minpoly)); e = 4: adds the ranks of (A - a)^1 for each
    root a and of f(A) for quadratic factors, which separates the squarefree classes used here"""
    e = len(A)
    if e == 2:
        # 2 x 2: the class is determined by the characteristic polynomial and whether A is scalar
        cp = tuple(charpoly(A, p))
        scalar = A[0][1] == 0 and A[1][0] == 0 and A[0][0] == A[1][1]
        return (cp, (A[0][0] * (p - 1) % p, 1) if scalar else cp)
    cp, mp = minpoly_key(A, p)
    if e <= 3:
        return (cp, mp)
    extra = []
    for a in range(p):
        extra.append(rank_mod_p(poly_eval_matrix([(-a) % p, 1], A, p), p))
    return (cp, mp, tuple(extra))


def companion(poly, p):
    """companion matrix of a monic polynomial (low-first, last coeff 1)"""
    d = len(poly) - 1
    C = [[0] * d for _ in range(d)]
    for i in range(1, d):
        C[i][i - 1] = 1
    for i in range(d):
        C[i][d - 1] = (-poly[i]) % p
    return C


def block_diag(*Ms):
    n = sum(len(M) for M in Ms)
    R = [[0] * n for _ in range(n)]
    o = 0
    for M in Ms:
        for i in range(len(M)):
            for j in range(len(M)):
                R[o + i][o + j] = M[i][j]
        o += len(M)
    return R


def gl_classes(e, p):
    """one representative per conjugacy class of GL(e, p) (e <= 3 complete; e = 4: squarefree charpoly only)"""
    reps = {}
    if e == 1:
        for a in range(1, p):
            reps[gl_key([[a]], p)] = [[a]]
    elif e == 2:
        for a in range(1, p):
            reps.setdefault(gl_key([[a, 0], [0, a]], p), [[a, 0], [0, a]])
            reps.setdefault(gl_key([[a, 1], [0, a]], p), [[a, 1], [0, a]])
            for b in range(a + 1, p):
                reps.setdefault(gl_key([[a, 0], [0, b]], p), [[a, 0], [0, b]])
        for f in monic_polys(2, p):
            reps.setdefault(gl_key(companion(f, p), p), companion(f, p))
    elif e == 3:
        for f in monic_polys(3, p):                      # cyclic matrices
            reps.setdefault(gl_key(companion(f, p), p), companion(f, p))
        for a in range(1, p):                             # invariant factors (x-a) | (x-a)(x-b): diag(a,a,b) and J2(a)+J1(a)
            for g in monic_polys(1, p):
                f2 = poly_mul([(-a) % p, 1], g, p)
                Mx = block_diag([[a]], companion(f2, p))
                reps.setdefault(gl_key(Mx, p), Mx)
            Mx = [[a if i == j else 0 for j in range(3)] for i in range(3)]
            reps.setdefault(gl_key(Mx, p), Mx)
    elif e == 4:
        for f in monic_polys(4, p):
            # squarefree test: gcd(f, f') = 1 via checking no repeated irreducible factor: use the derivative
            df = [(i * f[i]) % p for i in range(1, len(f))]
            if not any(df):
                continue
            if poly_gcd_is_one(f, df, p):
                reps.setdefault(gl_key(companion(f, p), p), companion(f, p))
    return list(reps.values())


def poly_gcd_is_one(a, b, p):
    def trim(x):
        while x and x[-1] == 0:
            x = x[:-1]
        return x
    a, b = trim(a[:]), trim(b[:])
    while b:
        # a mod b
        r = a[:]
        while len(r) >= len(b):
            f = r[-1] * pow(b[-1], -1, p) % p
            sh = len(r) - len(b)
            for i, c in enumerate(b):
                r[sh + i] = (r[sh + i] - f * c) % p
            r = trim(r[:-1]) if r and r[-1] == 0 else trim(r)
            if not r:
                break
        a, b = b, r
    return len(a) == 1


def gl2_all(p):
    for a, b, c, d in itertools.product(range(p), repeat=4):
        if (a * d - b * c) % p:
            yield [[a, b], [c, d]]


def eig1_fixed(A, p):
    """dimension of ker(A - I)"""
    e = len(A)
    Mx = [[(A[i][j] - (i == j)) % p for j in range(e)] for i in range(e)]
    return e - rank_mod_p(Mx, p)


# ---------------- block automorphism class lists ----------------
# each entry: dict(spec=<block spec>, order=o, key=<exact class key or None if twisted>, fix_pow=callable(m)->#fixed points
# of alpha^m, zfix_pow=callable(m)->#fixed points of alpha^m on Z(block)), built lazily and cached
_cache = {}


def perm_order(perm):
    n = len(perm); o = 1; cur = perm.copy(); idx = np.arange(n)
    while not np.array_equal(cur, idx):
        cur = perm[cur]; o += 1
    return o


def perm_fix_counts(perm, M):
    """[#fixed points of perm^m for m in 0..M-1] (perm^M = id assumed)"""
    n = len(perm); idx = np.arange(n); out = []
    cur = idx.copy()
    for m in range(M):
        out.append(int((cur == idx).sum()))
        cur = perm[cur]
    return out


def centre_mask(T):
    n = len(T)
    return np.all(T == T.T, axis=1)


def entry_from_perm(spec, perm, T, key, extra=None):
    o = perm_order(perm)
    zm = centre_mask(T)
    fx = perm_fix_counts(perm, o)
    zfx = []
    cur = np.arange(len(perm))
    for m in range(o):
        zfx.append(int(((cur == np.arange(len(perm))) & zm).sum())); cur = perm[cur]
    return {"spec": spec, "order": o, "key": key, "fix": fx, "zfix": zfx, "size": len(perm), "extra": extra or {}}


def autos_cyc(q, target=None):
    """units u mod q with u^target = 1 (all units when target is None): order, fixed points of u^m"""
    k = ("cyc", q, target)
    if k in _cache:
        return _cache[k]
    out = []
    for u in range(1, q):
        if gcd(u, q) != 1 or (target is not None and pow(u, target, q) != 1 % q):
            continue
        o = mult_order(u, q)
        fx = [gcd(pow(u, m, q) - 1, q) if q > 1 else 1 for m in range(o)]
        out.append({"spec": {"t": "cyc", "q": q, "u": u}, "order": o, "key": ("cyc", q, u), "fix": fx, "zfix": fx, "size": q, "extra": {}})
    _cache[k] = out
    return out


def autos_vec(p, e):
    k = ("vec", p, e)
    if k in _cache:
        return _cache[k]
    out = []
    for A in gl_classes(e, p):
        n_, T, perm, _ = grp.block_vec(p, A)
        out.append(entry_from_perm({"t": "vec", "p": p, "A": A}, perm.astype(np.int64), T, ("vec", p, e, gl_key(A, p))))
    _cache[k] = out
    return out


def mat_inv2(N, p):
    d = pow((N[0][0] * N[1][1] - N[0][1] * N[1][0]) % p, -1, p)
    return [[N[1][1] * d % p, (-N[0][1]) * d % p], [(-N[1][0]) * d % p, N[0][0] * d % p]]


def frob_B(p, k):
    """one fixed-point-free B of order k >= 3 (odd) in GL(2, p): the companion matrix of an irreducible quadratic whose
    root has order k when k | p^2 - 1 and k does not divide p - 1 (non-split torus); else diag(b1, b2) with b1 != b2,
    both != 1, lcm of orders = k (split); else None (scalar B is handled as the 'frobS' option)."""
    if k % 2 == 0 or k < 3 or (p * p - 1) % k:
        return None
    if (p - 1) % k:
        for f in monic_polys(2, p):
            B = companion(f, p)
            if grp.matorder(B, p) == k and eig1_fixed(B, p) == 0 and eig1_fixed(grp.matpow(B, 2, p), p) == 0:
                # irreducible iff no eigenvalue in F_p: check no root of f
                if all((a * a + f[1] * a + f[0]) % p for a in range(p)):
                    return B
        return None
    for b1 in range(2, p):
        for b2 in range(b1 + 1, p):
            if lcm(mult_order(b1, p), mult_order(b2, p)) == k:
                return [[b1, 0], [0, b2]]
    return None


def frob_normalizer_list(p, B, k):
    """all (N, j) with N B N^-1 = B^j, from the torus structure (no GL(2,p) scan): non-split: T = {aI + bB} and T F with
    F the Frobenius (F B F^-1 = B^p); split: monomial matrices"""
    Bp = [grp.matpow(B, i, p) for i in range(k)]
    out = []
    if B[0][1] == 0 and B[1][0] == 0:            # split: diagonal and anti-diagonal matrices
        cands = [[[a, 0], [0, d]] for a in range(1, p) for d in range(1, p)] + [[[0, b], [c, 0]] for b in range(1, p) for c in range(1, p)]
    else:
        Bpp = grp.matpow(B, p, p)                 # B^p = c0 I + c1 B  (B[1][0] == 1 for a companion matrix)
        c1 = Bpp[1][0] * pow(B[1][0], -1, p) % p
        c0 = (Bpp[0][0] - c1 * B[0][0]) % p
        F = [[1, c0], [0, c1]]
        T = [[[(a + b * B[0][0]) % p, (b * B[0][1]) % p], [(b * B[1][0]) % p, (a + b * B[1][1]) % p]] for a in range(p) for b in range(p) if (a, b) != (0, 0)]
        cands = T + [grp.matmul(t, F, p) for t in T]
    for N in cands:
        if (N[0][0] * N[1][1] - N[0][1] * N[1][0]) % p == 0:
            continue
        conj = grp.matmul(grp.matmul(N, B, p), mat_inv2(N, p), p)
        j = next((i for i in range(k) if Bp[i] == conj), None)
        if j is not None:
            out.append((N, j))
    return out


def autos_frob(p, k):
    ck = ("frob", p, k)
    if ck in _cache:
        return _cache[ck]
    B = frob_B(p, k)
    out = []
    if B is None:
        _cache[ck] = out
        return out
    Nrm = frob_normalizer_list(p, B, k)
    Ninvs = [mat_inv2(N, p) for N, _ in Nrm]
    Bp = [grp.matpow(B, i, p) for i in range(k)]
    # (A, e) up to conjugacy by the normalizer: canonical = min over N of flattened N A N^-1
    seen = {}
    for A, e in Nrm:
        canon = min(tuple(v for row in grp.matmul(grp.matmul(N, A, p), Ni, p) for v in row) for (N, _), Ni in zip(Nrm, Ninvs))
        seen.setdefault((e, canon), A)
    for (e, canon), A in seen.items():
        n_, T, perm, _ = grp.block_frob(p, B, A, e)
        out.append(entry_from_perm({"t": "frob", "p": p, "B": B, "A": A, "e": e}, perm.astype(np.int64), T, ("frob", p, k, e, canon)))
        # translation twists when A has eigenvalue 1: alpha' = alpha o inn((w, 0)), w over coset reps of (A - I) F_p^2
        if eig1_fixed(A, p) > 0 and A != [[1, 0], [0, 1]]:   # A = I twists are inner (central element in G)
            img = set()
            for w in itertools.product(range(p), repeat=2):
                v = ((A[0][0] - 1) * w[0] + A[0][1] * w[1]) % p, (A[1][0] * w[0] + (A[1][1] - 1) * w[1]) % p
                img.add(v)
            cosets = {}
            for w in itertools.product(range(p), repeat=2):
                key = min((((w[0] + v[0]) % p), ((w[1] + v[1]) % p)) for v in img)
                cosets.setdefault(key, w)
            for key, w in cosets.items():
                if key == (0, 0):
                    continue
                # inn((w,0)) as a permutation on the block: (v, i) -> (v + (I - B^i) w, i)
                n2 = p * p * k
                idx = np.arange(n2); i_ = idx % k; t = idx // k; b_ = t % p; a_ = t // p
                wa = np.zeros(n2, dtype=np.int64); wb = np.zeros(n2, dtype=np.int64)
                for i in range(k):
                    Bi = Bp[i]
                    m = i_ == i
                    wa[m] = (w[0] - (Bi[0][0] * w[0] + Bi[0][1] * w[1])) % p
                    wb[m] = (w[1] - (Bi[1][0] * w[0] + Bi[1][1] * w[1])) % p
                inn = (((a_ + wa) % p) * p + (b_ + wb) % p) * k + i_
                perm2 = perm.astype(np.int64)[inn]          # alpha o inn
                out.append(entry_from_perm({"t": "frob", "p": p, "B": B, "A": A, "e": e, "w": list(w)}, perm2, T, None, {"twist": True}))
    _cache[ck] = out
    return out


def autos_heis(p):
    ck = ("heis", p)
    if ck in _cache:
        return _cache[ck]
    out = []
    for A in gl_classes(2, p):
        n_, T, perm, _ = grp.block_heis(p, A)
        out.append(entry_from_perm({"t": "heis", "p": p, "A": A}, perm.astype(np.int64), T, ("heis", p, gl_key(A, p))))
        if eig1_fixed(A, p) > 0 and A != [[1, 0], [0, 1]]:
            # twists by inn(h), h over coset representatives of (A - I) F_p^2 in F_p^2 (inner automorphisms = conjugation by
            # (w, 0): (v, z) -> (v, z + w(w, v)))
            img = set()
            for w in itertools.product(range(p), repeat=2):
                img.add((((A[0][0] - 1) * w[0] + A[0][1] * w[1]) % p, (A[1][0] * w[0] + (A[1][1] - 1) * w[1]) % p))
            cosets = {}
            for w in itertools.product(range(p), repeat=2):
                key = min((((w[0] + v[0]) % p), ((w[1] + v[1]) % p)) for v in img)
                cosets.setdefault(key, w)
            n2 = p ** 3
            idx = np.arange(n2); z_ = idx % p; t = idx // p; b_ = t % p; a_ = t // p
            for key, w in cosets.items():
                if key == (0, 0):
                    continue
                om = (w[0] * b_ - w[1] * a_) % p
                inn = (a_ * p + b_) * p + (z_ + om) % p
                perm2 = perm.astype(np.int64)[inn]
                out.append(entry_from_perm({"t": "heis", "p": p, "A": A, "w": list(w)}, perm2, T, None, {"twist": True}))
    _cache[ck] = out
    return out


def autos_ab2(p, a, b):
    ck = ("ab2", p, a, b)
    if ck in _cache:
        return _cache[ck]
    qa, qb = p ** a, p ** b
    n_ = qa * qb
    idx = np.arange(n_); y_ = idx % qb; x_ = idx // qb
    # all automorphisms: images of e1 = (1,0) (order qa) and e2 = (0,1) (order qb) that generate
    autos = []
    for x1 in range(qa):
        for y1 in range(qb):
            for x2 in range(qa):
                for y2 in range(qb):
                    img = ((x1 * x_ + x2 * y_) % qa) * qb + (y1 * x_ + y2 * y_) % qb
                    if len(set(img.tolist())) != n_:
                        continue
                    # must be a homomorphism: automatically (linear) if well defined: need ord(img e1) | qa, ord(img e2) | qb
                    if (qa * x1) % qa or (qa * y1) % qb or (qb * x2) % qa or (qb * y2) % qb:
                        continue
                    autos.append((((x1, y1), (x2, y2)), img))
    # conjugacy classes (brute force on permutations)
    perms = [im for _, im in autos]
    invs = [np.argsort(im) for im in perms]
    classes = {}
    done = set()
    for i, (spec, im) in enumerate(autos):
        if i in done:
            continue
        cl = set()
        for j in range(len(perms)):
            conj = perms[j][im[invs[j]]]
            cl.add(tuple(conj.tolist()))
        for j2, (_, im2) in enumerate(autos):
            if tuple(im2.tolist()) in cl:
                done.add(j2)
        classes[i] = spec
    T = grp.block_ab2(p, a, b, [[1, 0], [0, 1]])[1]
    out = []
    for i, spec in classes.items():
        # no power-map reduction for ab2 blocks (the brute-force conjugacy test of alpha^j is too slow): kept like twisted entries
        out.append(entry_from_perm({"t": "ab2", "p": p, "a": a, "b": b, "img": [list(spec[0]), list(spec[1])]}, perms[i].astype(np.int64), T, ("ab2", p, a, b, i), {"twist": True}))
    _cache[ck] = out
    return out


def p_part_options(p, e):
    """block option lists for the p-part p^e (excluding frob blocks, which are handled at the K level)"""
    opts = [[("cyc", p ** e)]]
    if e == 2 and p <= 31:                           # p = 37, 59 occur only with M <= 15 (central Z_3 / inverted F_p^2): skipped
        opts.append([("vec", p, 2)])
    if e == 3:
        if p <= 7:                                   # ab2 by brute force over Aut: |Aut(Z_49 x Z_7)| = 12348 is the limit
            opts.append([("ab2", p, 2, 1)])
            opts.append([("vec", p, 3)]); opts.append([("heis", p)])
    if e == 4 and p == 3:                            # 81: Z_27 x Z_3, Z_9 x Z_9, F_3^4 (squarefree classes); p >= 5: cyclic only
        opts.append([("ab2", p, 3, 1)]); opts.append([("ab2", p, 2, 2)])
        opts.append([("vec", p, 4)])
    return opts


def get_autos(bl):
    t = bl[0]
    if t == "cyc":
        return autos_cyc(bl[1])
    if t == "vec":
        return autos_vec(bl[1], bl[2])
    if t == "heis":
        return autos_heis(bl[1])
    if t == "ab2":
        return autos_ab2(bl[1], bl[2], bl[3])
    if t == "frob":
        return autos_frob(bl[1], bl[2])
    raise ValueError(t)


def decompositions(K):
    """block lists covering |K| = K: standard p-part options, plus frob blocks F_p^2 x| C_k"""
    f = factor(K)
    results = []

    def rec(primes, acc):
        if not primes:
            results.append(list(acc)); return
        p = primes[0]; e = f[p]
        for opt in p_part_options(p, e):
            rec(primes[1:], acc + opt)
    rec(sorted(f), [])
    # frob blocks: p^2 | K, k | K/p^2, k odd >= 3, k | p^2 - 1, B fixed-point-free exists
    for p, e in f.items():
        if e < 2 or p > 17:          # p >= 19: no faithful C_M (M >= 4) in the window (normalizer orders 2(p^2-1) / 2(p-1)^2 vs n/(p^2 k))
            continue
        rest0 = K // (p * p)
        for k in range(3, p * p, 2):
            if (p * p - 1) % k or rest0 % k:
                continue
            if frob_B(p, k) is None:
                continue
            rest = rest0 // k
            fr = factor(rest)
            sub = []

            def rec2(primes, acc):
                if not primes:
                    sub.append(list(acc)); return
                pp = primes[0]; ee = fr[pp]
                for opt in p_part_options(pp, ee):
                    rec2(primes[1:], acc + opt)
            rec2(sorted(fr), [])
            for s in sub:
                results.append([("frob", p, k)] + s)
    return results


def main():
    LO, HI = int(sys.argv[1]), int(sys.argv[2])
    arg = lambda k, d: type(d)(sys.argv[sys.argv.index(k) + 1]) if k in sys.argv else d
    mmin, mmax, imax, zmax = arg("--mmin", 4), arg("--mmax", 150), arg("--imax", 31), arg("--zmax", 2)
    outfn = arg("--out", "")
    count_only = "--count-only" in sys.argv
    cand = []
    stats = defaultdict(int)
    import time
    t0 = time.time()
    for n in range(LO + 1, HI + 1):
        if n % 50 == 0:
            print(f"# n = {n}: {len(cand)} candidates so far, {time.time() - t0:.0f}s", file=sys.stderr, flush=True)
        for M in range(mmin, mmax + 1):
            if n % M:
                continue
            K = n // M
            if K % 2 == 0 or K < 3:
                continue
            for dec in decompositions(K):
                for variant in ("faithful", "kernel2"):
                    target = M if variant == "faithful" else (M // 2 if M % 2 == 0 else 0)
                    if target == 0:
                        continue
                    # block lists filtered by order | target (cyclic blocks are built per target: only units with u^target = 1)
                    fl = [[a for a in (autos_cyc(bl[1], target) if bl[0] == "cyc" else get_autos(bl)) if target % a["order"] == 0] for bl in dec]
                    if any(len(L) == 0 for L in fl):
                        continue
                    for combo in itertools.product(*fl):
                        o = 1
                        for a in combo:
                            o = lcm(o, a["order"])
                        if o != target:
                            continue
                        stats["combos"] += 1
                        # involutions: |K| / prod fix(alpha_i^(M/2)) for M even (faithful); kernel2: y^(M/2) central -> 1 involution
                        if M % 2 == 0:
                            if variant == "faithful":
                                fixp = 1
                                for a in combo:
                                    fixp *= a["fix"][(M // 2) % a["order"]]
                                ninv = K // fixp
                            else:
                                ninv = 1
                        else:
                            ninv = 0
                        if ninv > imax:
                            stats["inv_cut"] += 1; continue
                        # centre pre-filter: fixed points of alpha on Z(K) must be trivial (|Z(G)| <= zmax = 2 allows only the kernel involution)
                        zf = 1
                        for a in combo:
                            zf *= a["zfix"][1 % a["order"]] if a["order"] > 1 else a["zfix"][0]
                        if zf > 1:
                            stats["zfix_cut"] += 1; continue
                        cand.append({"n": n, "M": M, "K": K, "variant": variant, "ninv": ninv, "dec": dec,
                                     "blocks": [a["spec"] for a in combo], "keys": [a["key"] for a in combo], "orders": [a["order"] for a in combo],
                                     "twist": any(a["extra"].get("twist") for a in combo)})
    print(f"# combos {stats['combos']}, involution cut {stats['inv_cut']}, Z(K)-fixed cut {stats['zfix_cut']}, surviving {len(cand)}", file=sys.stderr)
    # power-map reduction: alpha ~ alpha^j for units j mod M, exact when no block is twisted: use key tuples of alpha_i^j,
    # computed by rebuilding each block's permutation power and matching its key among the block's class list
    kept = []
    seen = set()
    for c in cand:
        if c["twist"]:
            kept.append(c); continue
        M = c["M"]
        canon_forms = []
        for j in range(1, M):
            if gcd(j, M) != 1:
                continue
            canon_forms.append(tuple(pow_key(bl, sp, o, j) for bl, sp, o in zip(c["dec"], c["blocks"], c["orders"])))
        cf = (c["n"], c["M"], c["variant"], tuple(c["dec"]), min(canon_forms))
        if cf in seen:
            continue
        seen.add(cf)
        kept.append(c)
    print(f"# after power-map reduction: {len(kept)}", file=sys.stderr)
    if count_only:
        return
    if outfn:
        with open(outfn, "w") as f:
            for i, c in enumerate(kept):
                spec = {"M": c["M"], "blocks": c["blocks"]}
                f.write(json.dumps({"id": f"g_{c['n']}_{i}", "n": c["n"], "M": c["M"], "K": c["K"], "variant": c["variant"], "ninv_pred": c["ninv"],
                                    "dec": [list(d) for d in c["dec"]], "twist": c["twist"], "spec": spec}) + "\n")


_pk_cache = {}


def frob_normalizer(p, k):
    ck = ("frobN", p, k)
    if ck not in _pk_cache:
        B = frob_B(p, k)
        _pk_cache[ck] = [(N, mat_inv2(N, p)) for N, _ in frob_normalizer_list(p, B, k)]
    return _pk_cache[ck]


def pow_key(bl, spec, order, j):
    """exact conjugacy-class key of alpha_i^j for an untwisted block automorphism given by its spec (cached)"""
    t = bl[0]
    jj = j % order
    ck = (t, json.dumps(spec, sort_keys=True), jj)
    if ck in _pk_cache:
        return _pk_cache[ck]
    if t == "cyc":
        q = spec["q"]
        key = ("cyc", q, pow(spec["u"], jj, q) if q > 1 else 1)
    elif t == "vec":
        p = spec["p"]; A = grp.matpow(spec["A"], jj, p)
        key = ("vec", p, len(A), gl_key(A, p))
    elif t == "heis":
        p = spec["p"]; A = grp.matpow(spec["A"], jj, p)
        key = ("heis", p, gl_key(A, p))
    elif t == "frob":
        p, k = bl[1], bl[2]
        A = grp.matpow(spec["A"], jj, p); e = pow(spec.get("e", 1), jj, k) if k > 1 else 0
        Nrm = frob_normalizer(p, k)
        canon = min(tuple(v for row in grp.matmul(grp.matmul(N, A, p), Ninv, p) for v in row) for N, Ninv in Nrm)
        key = ("frob", p, k, e, canon)
    elif t == "ab2":
        p, a, b = bl[1], bl[2], bl[3]
        n_, T, perm, _ = grp.block_ab2(p, a, b, spec["img"])
        perm = perm.astype(np.int64)
        pj = np.arange(n_)
        for _ in range(jj):
            pj = perm[pj]
        key = None
        for ent in autos_ab2(p, a, b):
            n2, T2, perm2, _ = grp.build_block(ent["spec"])
            if _ab2_conj(bl, perm2.astype(np.int64), pj):
                key = ent["key"]; break
        if key is None:
            raise RuntimeError("ab2 class not found")
    else:
        raise ValueError(t)
    _pk_cache[ck] = key
    return key


_ab2_autos = {}


def _ab2_conj(bl, perm, pj):
    """is pj conjugate to perm under Aut(Z_{p^a} x Z_{p^b})?"""
    p, a, b = bl[1], bl[2], bl[3]
    ck = (p, a, b)
    if ck not in _ab2_autos:
        qa, qb = p ** a, p ** b
        n_ = qa * qb
        idx = np.arange(n_); y_ = idx % qb; x_ = idx // qb
        autos = []
        for x1 in range(qa):
            for y1 in range(qb):
                for x2 in range(qa):
                    for y2 in range(qb):
                        if (qa * x1) % qa or (qa * y1) % qb or (qb * x2) % qa or (qb * y2) % qb:
                            continue
                        img = ((x1 * x_ + x2 * y_) % qa) * qb + (y1 * x_ + y2 * y_) % qb
                        if len(set(img.tolist())) == n_:
                            autos.append(img)
        _ab2_autos[ck] = autos
    for g in _ab2_autos[ck]:
        gi = np.argsort(g)
        if np.array_equal(g[perm[gi]], pj):
            return True
    return False


if __name__ == "__main__":
    main()
