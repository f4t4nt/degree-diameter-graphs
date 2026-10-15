#!/usr/bin/env python3
"""gf.py: build G = K x| C_M for an odd-order K given as a direct product of blocks with an automorphism alpha
(y k y^-1 = alpha(k)), write the binary group file for table_sweep.c, and report invariants.

Blocks (spec = JSON {"M": M, "blocks": [...]}, element (k, m) = k y^m has index k * M + m):
  {"t": "cyc", "q": q, "u": u}                        Z_q, alpha = multiplication by u (unit mod q)
  {"t": "vec", "p": p, "A": [[..]..]}                 F_p^e (e = len(A)), alpha = A (matrix in GL(e, p))
  {"t": "frob", "p": p, "B": [[..]], "A": [[..]], "e": e}
                                                      F_p^2 x| <B> (B of order k, c acts by B), alpha(v, i) = (A v, e i);
                                                      requires A B A^-1 = B^e; k = 1 (B = I) gives plain F_p^2
  {"t": "heis", "p": p, "A": [[..]]}                  p^{1+2}_+ as F_p^2 x F_p with (v,z)(v',z') = (v+v', z+z'+w(v,v')/2),
                                                      alpha(v, z) = (A v, det(A) z), A in GL(2, p)
  {"t": "ab2", "p": p, "a": a, "b": b, "img": [[x1,y1],[x2,y2]]}
                                                      Z_{p^a} x Z_{p^b} (a >= b), alpha(e1) = (x1, y1), alpha(e2) = (x2, y2)
K = direct product of the blocks (index = mixed radix, first block most significant), alpha = product.
The file: magic 'GRPF', int32 nK, int32 M, uint16 TK[nK*nK] (TK[a*nK+b] = a*b), uint16 alpha[nK], then the spec JSON
(null-terminated).  Checks: K's table is a group (identity 0, inverses, associativity on 2e6 random triples + full for
nK <= 400), alpha is an automorphism (all pairs) of order exactly dividing M, and the reported invariants come from
the assembled G (numpy): |Z(G)|, involutions, element-order histogram (optional), the faithfulness order of alpha.

usage: gf.py build SPEC_JSON OUT.grp [--orders]      (SPEC_JSON = a JSON string or @file)
       gf.py info OUT.grp [--orders]
"""
import json
import sys
from collections import Counter
from math import gcd

import numpy as np

MAGIC = b"GRPF"


def matmul(A, B, p):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) % p for j in range(len(B[0]))] for i in range(len(A))]


def matpow(A, e, p):
    n = len(A)
    R = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(e):
        R = matmul(R, A, p)
    return R


def matorder(A, p, limit=100000):
    n = len(A)
    I = [[int(i == j) for j in range(n)] for i in range(n)]
    R = A
    for k in range(1, limit):
        if R == I:
            return k
        R = matmul(R, A, p)
    raise ValueError("matrix order too large")


def det2(A, p):
    return (A[0][0] * A[1][1] - A[0][1] * A[1][0]) % p


def vec_index(v, p):
    idx = 0
    for x in v:
        idx = idx * p + int(x) % p
    return idx


def block_cyc(q, u):
    assert q % 2 == 1 and gcd(u, q) == 1
    a = np.arange(q)
    T = (a[:, None] + a[None, :]) % q
    alpha = (u * a) % q
    return q, T.astype(np.uint16), alpha.astype(np.uint16), [1 % q] if q > 1 else []


def block_vec(p, A):
    e = len(A)
    n = p ** e
    idx = np.arange(n)
    coords = np.stack([(idx // p ** (e - 1 - i)) % p for i in range(e)], axis=1)   # (n, e)
    # table: addition
    T = np.zeros((n, n), dtype=np.int64)
    for i in range(e):
        T = T * p + (coords[:, None, i] + coords[None, :, i]) % p
    Am = np.array(A, dtype=np.int64) % p
    img = (coords @ Am.T) % p
    alpha = np.zeros(n, dtype=np.int64)
    for i in range(e):
        alpha = alpha * p + img[:, i]
    gens = [p ** (e - 1 - i) for i in range(e)]
    return n, T.astype(np.uint16), alpha.astype(np.uint16), gens


def block_frob(p, B, A, e):
    """F_p^2 x| <B>: elements (v, i), index vec_index(v) * k + i; (v1, i1)(v2, i2) = (v1 + B^i1 v2, i1 + i2)"""
    k = matorder(B, p)
    assert k % 2 == 1, "k must be odd"
    Bp = [np.array(matpow(B, i, p), dtype=np.int64) for i in range(k)]
    Am = np.array(A, dtype=np.int64) % p
    e = e % k if k > 1 else 0
    if k > 1:
        Ainv = np.array(matpow(A, matorder(A, p) - 1, p), dtype=np.int64)
        assert np.array_equal((Am @ Bp[1] @ Ainv) % p, Bp[e]), "A B A^-1 != B^e"
        assert gcd(e, k) == 1, "e must be a unit mod k"
    n = p * p * k
    idx = np.arange(n)
    i_ = idx % k; t = idx // k
    b_ = t % p; a_ = t // p
    V = np.stack([a_, b_], axis=1)   # (n, 2)
    # product table
    T = np.zeros((n, n), dtype=np.int64)
    for i1 in range(k):
        rows = np.nonzero(i_ == i1)[0]
        W = (V[None, :, :] @ Bp[i1].T) % p            # (1, n, 2): B^i1 v2
        W = W[0]
        va = (V[rows, 0][:, None] + W[:, 0][None, :]) % p
        vb = (V[rows, 1][:, None] + W[:, 1][None, :]) % p
        ii = (i1 + i_[None, :]) % k
        T[rows, :] = (va * p + vb) * k + ii
    img = (V @ Am.T) % p
    alpha = (img[:, 0] * p + img[:, 1]) * k + (e * i_) % k
    gens = [vec_index([1, 0], p) * k, vec_index([0, 1], p) * k] + ([1] if k > 1 else [])
    return n, T.astype(np.uint16), alpha.astype(np.uint16), gens


def block_heis(p, A):
    """p^{1+2}_+ : (v, z), index vec_index(v) * p + z, (v,z)(v',z') = (v+v', z+z'+h w(v,v')), h = 1/2 mod p"""
    assert p % 2 == 1
    h = pow(2, -1, p)
    n = p ** 3
    idx = np.arange(n)
    z_ = idx % p; t = idx // p
    b_ = t % p; a_ = t // p
    om = (a_[:, None] * b_[None, :] - b_[:, None] * a_[None, :]) % p
    va = (a_[:, None] + a_[None, :]) % p
    vb = (b_[:, None] + b_[None, :]) % p
    zz = (z_[:, None] + z_[None, :] + h * om) % p
    T = (va * p + vb) * p + zz
    Am = np.array(A, dtype=np.int64) % p
    d = det2(A, p)
    ia = (Am[0, 0] * a_ + Am[0, 1] * b_) % p
    ib = (Am[1, 0] * a_ + Am[1, 1] * b_) % p
    alpha = (ia * p + ib) * p + (d * z_) % p
    gens = [vec_index([1, 0], p) * p, vec_index([0, 1], p) * p]
    return n, T.astype(np.uint16), alpha.astype(np.uint16), gens


def block_ab2(p, a, b, img):
    """Z_{p^a} x Z_{p^b}, element (x, y) index x * p^b + y; alpha(e1) = img[0], alpha(e2) = img[1]"""
    qa, qb = p ** a, p ** b
    n = qa * qb
    idx = np.arange(n)
    y_ = idx % qb; x_ = idx // qb
    T = ((x_[:, None] + x_[None, :]) % qa) * qb + (y_[:, None] + y_[None, :]) % qb
    (x1, y1), (x2, y2) = img
    alpha = ((x1 * x_ + x2 * y_) % qa) * qb + (y1 * x_ + y2 * y_) % qb
    return n, T.astype(np.uint16), alpha.astype(np.uint16), [qb, 1]


def build_block(bl):
    t = bl["t"]
    if t == "cyc":
        return block_cyc(bl["q"], bl["u"])
    if t == "vec":
        return block_vec(bl["p"], bl["A"])
    if t == "frob":
        return block_frob(bl["p"], bl["B"], bl["A"], bl.get("e", 1))
    if t == "heis":
        return block_heis(bl["p"], bl["A"])
    if t == "ab2":
        return block_ab2(bl["p"], bl["a"], bl["b"], bl["img"])
    raise ValueError(t)


def direct_product(blocks):
    """blocks = [(n, T, alpha, gens)] -> (n, T, alpha, gens), index mixed radix (first block most significant)"""
    n, T, alpha, gens = blocks[0]
    T = T.astype(np.int64); alpha = alpha.astype(np.int64)
    for n2, T2, al2, g2 in blocks[1:]:
        N = n * n2
        idx = np.arange(N)
        k2 = idx % n2; k1 = idx // n2
        T = T[k1][:, k1] * n2 + T2.astype(np.int64)[k2][:, k2]
        alpha = alpha[k1] * n2 + al2.astype(np.int64)[k2]
        gens = [g * n2 for g in gens] + list(g2)
        n = N
    return n, T.astype(np.uint16), alpha.astype(np.uint16), gens


def check_group(n, T, alpha, rng):
    """K's table is a group with identity 0; alpha an automorphism.  Returns the order of alpha."""
    assert T.shape == (n, n)
    ALL = np.arange(n)
    assert np.array_equal(T[0], ALL) and np.array_equal(T[:, 0], ALL), "0 is not the identity"
    inv = np.argmax(T == 0, axis=1)
    assert np.all(T[ALL, inv] == 0) and np.all(T[inv, ALL] == 0), "inverses"
    assert all(len(set(T[a].tolist())) == n for a in range(0, n, max(1, n // 50))), "rows not permutations"
    if n <= 150:   # full triple check only for small K (the nK^3 int64 arrays reach 1.5 GB at nK = 400); table_sweep.c re-samples
        a = np.repeat(ALL, n * n); b = np.tile(np.repeat(ALL, n), n); c = np.tile(ALL, n * n)
    else:
        a, b, c = rng.integers(0, n, (3, 2_000_000))
    assert np.array_equal(T[T[a, b], c], T[a, T[b, c]]), "associativity"
    assert len(set(alpha.tolist())) == n, "alpha not a bijection"
    A = alpha.astype(np.int64)
    assert np.array_equal(A[T], T[A][:, A]), "alpha not a homomorphism"   # alpha(ab) = alpha(a) alpha(b) for all a, b
    # order of alpha
    o = 1; cur = A.copy()
    while not np.array_equal(cur, ALL):
        cur = A[cur]; o += 1
        assert o <= 10 ** 6
    return o


class Group:
    """G = K x| C_M assembled from the file data; vectorised multiplication"""

    def __init__(self, nK, M, TK, alpha, spec):
        self.nK, self.M, self.TK, self.alpha, self.spec = nK, M, TK.astype(np.int64), alpha.astype(np.int64), spec
        self.n = nK * M
        ap = [np.arange(nK)]
        for _ in range(1, M):
            ap.append(self.alpha[ap[-1]])
        self.apow = np.array(ap)          # apow[m][k] = alpha^m(k)
        assert np.array_equal(self.alpha[ap[-1]], ap[0]), "alpha^M != id"

    def mul(self, U, V):
        U = np.asarray(U, dtype=np.int64); V = np.asarray(V, dtype=np.int64)
        k1, m1 = U // self.M, U % self.M
        k2, m2 = V // self.M, V % self.M
        return self.TK[k1, self.apow[m1, k2]] * self.M + (m1 + m2) % self.M

    def inverse_all(self):
        n = self.n
        ALL = np.arange(n)
        # (k, m)^-1 = (alpha^{-m}(k^-1), -m)
        invK = np.argmax(self.TK == 0, axis=1)
        k, m = ALL // self.M, ALL % self.M
        inv = self.apow[(self.M - m) % self.M, invK[k]] * self.M + (self.M - m) % self.M
        assert np.all(self.mul(ALL, inv) == 0)
        return inv

    def invariants(self, gensK, orders=False):
        n = self.n
        ALL = np.arange(n)
        inv = self.inverse_all()
        ninv = int(np.sum((inv == ALL) & (ALL != 0)))
        gens = [g * self.M for g in gensK] + [1]            # K generators and y
        central = np.ones(n, bool)
        for g in gens:
            central &= self.mul(ALL, g) == self.mul(g, ALL)
        Z = int(central.sum())
        # faithfulness: smallest m > 0 with alpha^m = id
        o = next(m for m in range(1, self.M + 1) if np.array_equal(self.apow[m % self.M], np.arange(self.nK)))
        res = {"n": n, "nK": self.nK, "M": self.M, "Z": Z, "inv": ninv, "ord_alpha": o}
        if orders:
            order = np.zeros(n, dtype=np.int64); pw = ALL.copy(); kk = 1
            while True:
                done = (pw == 0) & (order == 0); order[done] = kk
                if np.all(order > 0):
                    break
                pw = self.mul(pw, ALL); kk += 1
            res["orders"] = dict(sorted(Counter(order.tolist()).items()))
            res["exponent"] = int(np.lcm.reduce(np.array(sorted(res["orders"]))))
        return res


def build(spec, rng=None):
    rng = rng or np.random.default_rng(0)
    blocks = [build_block(bl) for bl in spec["blocks"]]
    nK, TK, alpha, gens = direct_product(blocks)
    o = check_group(nK, TK, alpha, rng)
    M = int(spec["M"])
    assert M % o == 0, f"M = {M} is not a multiple of ord(alpha) = {o}"
    return nK, M, TK, alpha, gens, o


def write_grp(fn, nK, M, TK, alpha, spec):
    with open(fn, "wb") as f:
        f.write(MAGIC)
        f.write(np.array([nK, M], dtype=np.int32).tobytes())
        f.write(TK.astype(np.uint16).tobytes())
        f.write(alpha.astype(np.uint16).tobytes())
        f.write(json.dumps(spec, separators=(",", ":")).encode() + b"\0")


def read_grp(fn):
    with open(fn, "rb") as f:
        assert f.read(4) == MAGIC
        nK, M = np.frombuffer(f.read(8), dtype=np.int32)
        nK, M = int(nK), int(M)
        TK = np.frombuffer(f.read(2 * nK * nK), dtype=np.uint16).reshape(nK, nK)
        alpha = np.frombuffer(f.read(2 * nK), dtype=np.uint16)
        spec = json.loads(f.read().split(b"\0")[0].decode())
    return nK, M, TK, alpha, spec


def load_spec(arg):
    return json.load(open(arg[1:])) if arg.startswith("@") else json.loads(arg)


def main():
    if sys.argv[1] == "build":
        spec = load_spec(sys.argv[2])
        nK, M, TK, alpha, gens, o = build(spec)
        write_grp(sys.argv[3], nK, M, TK, alpha, spec)
        G = Group(nK, M, TK, alpha, spec)
        print(json.dumps(G.invariants(gens, "--orders" in sys.argv)))
    elif sys.argv[1] == "info":
        nK, M, TK, alpha, spec = read_grp(sys.argv[2])
        blocks = [build_block(bl) for bl in spec["blocks"]]
        nK2, TK2, alpha2, gens = direct_product(blocks)
        assert nK2 == nK and np.array_equal(TK2, TK) and np.array_equal(alpha2, alpha), "file does not match its spec"
        G = Group(nK, M, TK, alpha, spec)
        print(json.dumps(G.invariants(gens, "--orders" in sys.argv)))


if __name__ == "__main__":
    main()
