#!/usr/bin/env python3
"""tfgroup.py: TFG1 group files for tablefree_sweep.c (table-free evaluator of G = K x| C_M) from gf.py specs, plus a standalone
Python implementation of the block group law (no gf.py / numpy tables) used by tf_rebuild.py to rebuild edge lists
independently of tablefree_sweep.c.

TFG1 file (binary): magic 'TFG1', int32 M, int32 nblocks, then per block (spec order):
  int32 type (0 cyc, 1 vec, 2 frob, 3 heis, 4 ab2), int32 size, int32 ngens, int32 naut,
  cyc:   int32 q, int32 u  (the block's alpha is x -> u x; gens = [1])
  other: uint16 table[size*size] (T[a*size+b] = a*b), uint16 alpha[size], int32 gens[ngens], uint16 aut[naut*size]
         (aut = automorphisms of the block commuting with its alpha, from autK.block_autos; identity included)
then the JSON spec (null-terminated).  Element indices are gf.py's: k*M + m with k the mixed-radix block index
(first block most significant), so sets are interchangeable with table_sweep.c for groups that fit a table.
Twists ("w" in frob / heis specs) are applied exactly as enum_family.py builds them (alpha o inn((w, 0))); gf.py's
builder ignores them, so twisted specs built by gf.py are NOT the same group as the TFG1 file.

usage: tfgroup.py build SPEC_JSON|@file OUT.tfg [--autcap 4000] [--noaut]
       tfgroup.py check SPEC_JSON|@file            (standalone law vs gf.py's table on a sample; small groups)
"""
import itertools
import json
import sys

import numpy as np

import gf  # noqa: E402  (block builders)

TYPES = {"cyc": 0, "vec": 1, "frob": 2, "heis": 3, "ab2": 4}


# ---------------------------------------------------------------- standalone block law (tuples, no tables) --------------
def mat_vec(X, v, p):
    return tuple(sum(X[i][k] * v[k] for k in range(len(v))) % p for i in range(len(X)))


def mat_mul(X, Y, p):
    n = len(X)
    return tuple(tuple(sum(X[i][k] * Y[k][j] for k in range(n)) % p for j in range(n)) for i in range(n))


def mat_pow(X, e, p):
    n = len(X)
    R = tuple(tuple(int(i == j) for j in range(n)) for i in range(n))
    for _ in range(e):
        R = mat_mul(R, X, p)
    return R


def mat_order(X, p):
    n = len(X)
    I = tuple(tuple(int(i == j) for j in range(n)) for i in range(n))
    R = X; o = 1
    while R != I:
        R = mat_mul(R, X, p); o += 1
        assert o < 10 ** 6
    return o


class Block:
    """one block: elements are tuples; methods mul, inv, alpha, index/unindex (gf.py's block index), gens, identity"""

    def __init__(self, bl):
        self.t = bl["t"]
        self.bl = bl
        if self.t == "cyc":
            self.q, self.u = bl["q"], bl["u"]
            self.size = self.q
        elif self.t == "vec":
            self.p = bl["p"]; self.A = tuple(tuple(r) for r in bl["A"]); self.e = len(self.A)
            self.size = self.p ** self.e
        elif self.t == "frob":
            self.p = bl["p"]; self.B = tuple(tuple(r) for r in bl["B"]); self.A = tuple(tuple(r) for r in bl["A"])
            self.k = mat_order(self.B, self.p)
            self.ee = bl.get("e", 1) % self.k if self.k > 1 else 0
            self.Bp = [mat_pow(self.B, i, self.p) for i in range(self.k)]
            self.w = tuple(bl["w"]) if "w" in bl else None
            self.size = self.p * self.p * self.k
        elif self.t == "heis":
            self.p = bl["p"]; self.A = tuple(tuple(r) for r in bl["A"])
            self.h = pow(2, -1, self.p)
            self.d = (self.A[0][0] * self.A[1][1] - self.A[0][1] * self.A[1][0]) % self.p
            self.w = tuple(bl["w"]) if "w" in bl else None
            self.size = self.p ** 3
        elif self.t == "ab2":
            self.p, self.a, self.b = bl["p"], bl["a"], bl["b"]
            self.qa, self.qb = self.p ** self.a, self.p ** self.b
            (self.x1, self.y1), (self.x2, self.y2) = bl["img"]
            self.size = self.qa * self.qb
        else:
            raise ValueError(self.t)

    def identity(self):
        if self.t == "cyc":
            return 0
        if self.t == "vec":
            return (0,) * self.e
        if self.t == "frob":
            return ((0, 0), 0)
        if self.t == "heis":
            return ((0, 0), 0)
        return (0, 0)

    def mul(self, x, y):
        if self.t == "cyc":
            return (x + y) % self.q
        if self.t == "vec":
            return tuple((a + b) % self.p for a, b in zip(x, y))
        if self.t == "frob":
            (v1, i1), (v2, i2) = x, y
            w = mat_vec(self.Bp[i1], v2, self.p)
            return (((v1[0] + w[0]) % self.p, (v1[1] + w[1]) % self.p), (i1 + i2) % self.k)
        if self.t == "heis":
            (v1, z1), (v2, z2) = x, y
            om = (v1[0] * v2[1] - v1[1] * v2[0]) % self.p
            return (((v1[0] + v2[0]) % self.p, (v1[1] + v2[1]) % self.p), (z1 + z2 + self.h * om) % self.p)
        return ((x[0] + y[0]) % self.qa, (x[1] + y[1]) % self.qb)

    def inv(self, x):
        if self.t == "cyc":
            return (-x) % self.q
        if self.t == "vec":
            return tuple((-a) % self.p for a in x)
        if self.t == "frob":
            v, i = x
            j = (-i) % self.k
            w = mat_vec(self.Bp[j], v, self.p)      # (v, i)^-1 = (-B^-i v, -i)
            return (((-w[0]) % self.p, (-w[1]) % self.p), j)
        if self.t == "heis":
            v, z = x
            return (((-v[0]) % self.p, (-v[1]) % self.p), (-z) % self.p)
        return ((-x[0]) % self.qa, (-x[1]) % self.qb)

    def alpha(self, x):
        if self.t == "cyc":
            return (self.u * x) % self.q
        if self.t == "vec":
            return mat_vec(self.A, x, self.p)
        if self.t == "frob":
            if self.w is not None:          # alpha o inn((w, 0)): (v, i) -> (v + (I - B^i) w, i) first
                v, i = x
                Bw = mat_vec(self.Bp[i], self.w, self.p)
                x = (((v[0] + self.w[0] - Bw[0]) % self.p, (v[1] + self.w[1] - Bw[1]) % self.p), i)
            v, i = x
            return (mat_vec(self.A, v, self.p), (self.ee * i) % self.k)
        if self.t == "heis":
            if self.w is not None:          # inn((w, 0)): (v, z) -> (v, z + w(w, v)) with w(w, v) = w0 v1 - w1 v0
                v, z = x
                x = (v, (z + self.w[0] * v[1] - self.w[1] * v[0]) % self.p)
            v, z = x
            return (mat_vec(self.A, v, self.p), (self.d * z) % self.p)
        x0, y0 = x
        return ((self.x1 * x0 + self.x2 * y0) % self.qa, (self.y1 * x0 + self.y2 * y0) % self.qb)

    def index(self, x):
        if self.t == "cyc":
            return x
        if self.t == "vec":
            idx = 0
            for a in x:
                idx = idx * self.p + a
            return idx
        if self.t == "frob":
            v, i = x
            return (v[0] * self.p + v[1]) * self.k + i
        if self.t == "heis":
            v, z = x
            return (v[0] * self.p + v[1]) * self.p + z
        return x[0] * self.qb + x[1]

    def unindex(self, idx):
        if self.t == "cyc":
            return idx
        if self.t == "vec":
            digs = []
            for _ in range(self.e):
                digs.append(idx % self.p); idx //= self.p
            return tuple(reversed(digs))
        if self.t == "frob":
            i = idx % self.k; t = idx // self.k
            return ((t // self.p, t % self.p), i)
        if self.t == "heis":
            z = idx % self.p; t = idx // self.p
            return ((t // self.p, t % self.p), z)
        return (idx // self.qb, idx % self.qb)

    def gens(self):
        if self.t == "cyc":
            return [1 % self.q] if self.q > 1 else []
        if self.t == "vec":
            return [tuple(int(i == j) for j in range(self.e)) for i in range(self.e)]
        if self.t == "frob":
            return [((1, 0), 0), ((0, 1), 0)] + ([((0, 0), 1)] if self.k > 1 else [])
        if self.t == "heis":
            return [((1, 0), 0), ((0, 1), 0)]
        return [(1, 0), (0, 1)]


class Law:
    """G = K x| C_M with the standalone block law; elements (tuple of block elements, m); gf.py indices"""

    def __init__(self, spec):
        self.M = int(spec["M"])
        self.blocks = [Block(bl) for bl in spec["blocks"]]
        self.sizes = [b.size for b in self.blocks]
        self.nK = 1
        for s in self.sizes:
            self.nK *= s
        self.n = self.nK * self.M

    def alpha_pow(self, k, m):
        out = []
        for b, x in zip(self.blocks, k):
            for _ in range(m):
                x = b.alpha(x)
            out.append(x)
        return tuple(out)

    def mul(self, g, h):
        k1, m1 = g; k2, m2 = h
        k2 = self.alpha_pow(k2, m1)
        return (tuple(b.mul(x, y) for b, x, y in zip(self.blocks, k1, k2)), (m1 + m2) % self.M)

    def inv(self, g):
        k, m = g
        ki = tuple(b.inv(x) for b, x in zip(self.blocks, k))
        return (self.alpha_pow(ki, (self.M - m) % self.M), (self.M - m) % self.M)

    def identity(self):
        return (tuple(b.identity() for b in self.blocks), 0)

    def index(self, g):
        k, m = g
        idx = 0
        for b, x in zip(self.blocks, k):
            idx = idx * b.size + b.index(x)
        return idx * self.M + m

    def unindex(self, idx):
        m = idx % self.M; k = idx // self.M
        digs = []
        for b in reversed(self.blocks):
            digs.append(b.unindex(k % b.size)); k //= b.size
        return (tuple(reversed(digs)), m)

    def elements(self):
        for digs in itertools.product(*[range(s) for s in self.sizes]):
            k = tuple(b.unindex(d) for b, d in zip(self.blocks, digs))
            for m in range(self.M):
                yield (k, m)

    def show(self, g):
        k, m = g
        parts = []
        for b, x in zip(self.blocks, k):
            if b.t == "cyc":
                parts.append(f"x{b.q}^{x}")
            elif b.t == "vec":
                parts.append("v" + str(list(x)))
            elif b.t == "frob":
                parts.append(f"(v{list(x[0])},c^{x[1]})")
            elif b.t == "heis":
                parts.append(f"(v{list(x[0])},z^{x[1]})")
            else:
                parts.append(f"ab({x[0]},{x[1]})")
        return "(" + ",".join(parts) + f") y^{m}"


# ---------------------------------------------------------------- TFG1 writer --------------------------------------------
def twisted_alpha(bl, perm):
    """apply enum_family's translation twist to a block alpha permutation: alpha' = alpha o inn((w, 0))"""
    if "w" not in bl:
        return perm
    w = bl["w"]; p = bl["p"]
    perm = perm.astype(np.int64)
    if bl["t"] == "frob":
        k = gf.matorder(bl["B"], p)
        Bp = [gf.matpow(bl["B"], i, p) for i in range(k)]
        n2 = p * p * k
        idx = np.arange(n2); i_ = idx % k; t = idx // k; b_ = t % p; a_ = t // p
        wa = np.zeros(n2, dtype=np.int64); wb = np.zeros(n2, dtype=np.int64)
        for i in range(k):
            Bi = Bp[i]; m = i_ == i
            wa[m] = (w[0] - (Bi[0][0] * w[0] + Bi[0][1] * w[1])) % p
            wb[m] = (w[1] - (Bi[1][0] * w[0] + Bi[1][1] * w[1])) % p
        inn = (((a_ + wa) % p) * p + (b_ + wb) % p) * k + i_
        return perm[inn]
    if bl["t"] == "heis":
        n2 = p ** 3
        idx = np.arange(n2); z_ = idx % p; t = idx // p; b_ = t % p; a_ = t // p
        om = (w[0] * b_ - w[1] * a_) % p
        inn = (a_ * p + b_) * p + (z_ + om) % p
        return perm[inn]
    raise ValueError("twist on " + bl["t"])


def block_auts(bl, cap):
    """automorphisms of the block commuting with its (possibly twisted) alpha: autK.block_autos on the untwisted spec, then
    filtered numerically against the twisted alpha (twists are rare; the filter keeps a subgroup)"""
    import aut_k as autK
    bl0 = {k: v for k, v in bl.items() if k != "w"}
    perms, label = autK.block_autos(bl0, cap)
    if "w" in bl:
        nb, T, alpha, _ = gf.build_block(bl0)
        al = twisted_alpha(bl, alpha)
        perms = [pm for pm in perms if np.array_equal(pm[al], al[pm])]
        label += f",twist-filtered{len(perms)}"
    return perms, label


_block_cache = {}


def built_block(bl, bi, law, autcap, noaut):
    """(T, alpha, gens, auts, label) of a non-cyclic block, cached per block spec (blocks repeat across the groups of a window)"""
    key = (json.dumps(bl, sort_keys=True), autcap, noaut)
    if key in _block_cache:
        return _block_cache[key]
    nb, T, alpha, gens = gf.build_block({k: v for k, v in bl.items() if k != "w"})
    alpha = twisted_alpha(bl, alpha)
    # cross-check the block table and alpha against the standalone law on a sample (all pairs when small)
    B = law.blocks[bi]
    rng = np.random.default_rng(bi)
    pairs = [(a, b) for a in range(nb) for b in range(nb)] if nb <= 200 else [tuple(rng.integers(0, nb, 2)) for _ in range(3000)]
    for a, b in pairs:
        assert int(T[a, b]) == B.index(B.mul(B.unindex(int(a)), B.unindex(int(b)))), "standalone law != gf table"
    for a in range(nb):
        assert int(alpha[a]) == B.index(B.alpha(B.unindex(a))), "standalone alpha != gf alpha"
    if noaut:
        auts = [np.arange(nb)]; lab = "noaut"
    else:
        auts, lab = block_auts(bl, autcap)
        if len(auts) > autcap:
            auts = [np.arange(nb)]; lab += ",dropped>cap"
        elif not any(np.array_equal(pm, np.arange(nb)) for pm in auts):
            auts = [np.arange(nb)] + auts
    res = (nb, T, alpha, gens, auts, lab)
    if len(_block_cache) > 400:
        _block_cache.clear()
    _block_cache[key] = res
    return res


def build(spec, out, autcap=4000, noaut=False):
    M = int(spec["M"])
    blocks = spec["blocks"]
    law = Law(spec)
    with open(out, "wb") as f:
        f.write(b"TFG1")
        f.write(np.array([M, len(blocks)], dtype=np.int32).tobytes())
        labels = []
        for bi, bl in enumerate(blocks):
            t = TYPES[bl["t"]]
            if bl["t"] == "cyc":
                q, u = bl["q"], bl["u"] % bl["q"]
                f.write(np.array([t, q, 1, 0], dtype=np.int32).tobytes())
                f.write(np.array([q, u], dtype=np.int32).tobytes())
                labels.append(f"cyc{q}")
                continue
            nb, T, alpha, gens, auts, lab = built_block(bl, bi, law, autcap, noaut)
            labels.append(lab)
            f.write(np.array([t, nb, len(gens), len(auts)], dtype=np.int32).tobytes())
            f.write(np.asarray(T, dtype=np.uint16).tobytes())
            f.write(np.asarray(alpha, dtype=np.uint16).tobytes())
            f.write(np.array(gens, dtype=np.int32).tobytes())
            f.write(np.stack([np.asarray(pm, dtype=np.uint16) for pm in auts]).tobytes())
        f.write(json.dumps(spec, separators=(",", ":")).encode() + b"\0")
    return labels


def check(spec):
    """standalone law against gf.py's assembled group on a sample of products (groups small enough for gf.py)"""
    nK, M, TK, alpha, gens, o = gf.build(spec)
    G = gf.Group(nK, M, TK, alpha, spec)
    law = Law(spec)
    assert law.n == G.n
    rng = np.random.default_rng(1)
    U = rng.integers(0, G.n, 5000); V = rng.integers(0, G.n, 5000)
    P = G.mul(U, V)
    for u, v, pr in zip(U.tolist(), V.tolist(), P.tolist()):
        assert law.index(law.mul(law.unindex(u), law.unindex(v))) == pr, "law mismatch"
    for u in U.tolist()[:500]:
        assert law.mul(law.unindex(u), law.inv(law.unindex(u))) == law.identity()
    return G.n


def load_spec(arg):
    return json.load(open(arg[1:])) if arg.startswith("@") else json.loads(arg)


if __name__ == "__main__":
    if sys.argv[1] == "build":
        spec = load_spec(sys.argv[2])
        cap = int(sys.argv[sys.argv.index("--autcap") + 1]) if "--autcap" in sys.argv else 4000
        labels = build(spec, sys.argv[3], cap, "--noaut" in sys.argv)
        print(json.dumps({"out": sys.argv[3], "blocks": labels}))
    elif sys.argv[1] == "check":
        spec = load_spec(sys.argv[2])
        print("standalone law == gf.py on 5000 products, n =", check(spec))
