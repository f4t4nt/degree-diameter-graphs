#!/usr/bin/env python3
"""Every metacyclic group G(N, M, r, t) of the requested order parity = <x, y | x^N = 1, y^M = x^t, y x y^-1 = x^r> (r^M = 1 mod N,
t (r - 1) = 0 mod N) with |G| = N M in [lo, hi], one presentation per proven equivalence (rules: r minimal in its class {r^l : l unit mod M}; t one per orbit under t -> t u, t -> t + j S_M, t -> t l;
split presentations whose N has a central prime-power factor coprime to M dropped, the group reappears with the factor
inside Z_M).  Abelian groups (r = 1) are listed only with --all.

Lift bound (rigorous): a 2-vertex lift with loop voltages a, b and link voltages e, c over G covers the same lift over
G^ab = G / <x^(r-1)> (order M c, c = gcd(N, r - 1)); there a walk of length <= D from (0, e) ends at (v, a^i b^j c^l)
with |i| + |j| + |l| <= D, so the ball has at most 2 L3(D) vertices, L3(D) = lattice points of the l1-ball of Z^3
(575 for D = 7, 833 for D = 8), and diameter <= D forces 2 M c <= 2 L3(D).  For the 4-link (Haar) shape the ball is
at most 2 L3(D) as well (exponents of c1, c2, c3).  Groups failing the bound are listed with shape flag '-'.

Columns: n N M r t c Mc ord |A| fam_B2a fam_haar fp shapes, where
  |A|      = exact size of the reduction group A = {x -> x^i, y -> x^j y^l : i unit mod N, l unit mod M, r^l = r,
             j S_M = t (i - l) mod N} (automorphisms of G; all of Aut(G) for prime N and a faithful action),
  fam_B2a  = n^3 / (8 |A|), fam_haar = n^3 / (12 |A|): rough reduced-family sizes (gen_b2.py gives the exact counts),
  fp       = isomorphism fingerprint: |Z(G)|, |G'|, exponent, element-order histogram (equal fp = possibly isomorphic;
             different fp = certainly non-isomorphic).
usage: groups_meta.py lo hi D odd|even|all [--all]
Column ninv = number of involutions (from the element table); the semi-edge shapes need ninv >= 1 (S3, S1a), >= 2 (S2a,
S2b) or >= 3 (S1b, S1c) and are flagged accordingly (no bound applied beyond non-abelian).
"""
import math
import sys

import numpy as np


def factor(n):
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


def phi(n):
    r = n
    for p in factor(n):
        r = r // p * (p - 1)
    return r


def mult_order(r, n):
    o, p = 1, r % n
    while p != 1:
        p = p * r % n
        o += 1
    return o


def l1_ball3(D):
    return sum(math.comb(3, k) * math.comb(D, k) * 2**k for k in range(4))


_ug = {}


def unit_gens(N):
    if N in _ug:
        return _ug[N]
    H = {1 % N}
    gens = []
    for g in range(2, N):
        if math.gcd(g, N) != 1 or g in H:
            continue
        gens.append(g)
        stack = list(H)
        while stack:
            v = stack.pop()
            w = v * g % N
            if w not in H:
                H.add(w)
                stack.append(w)
    _ug[N] = gens
    return gens


def t_orbit_reps(N, M, r, c, SM, ls):
    step = N // c
    gens = unit_gens(N) + ls
    seen = set()
    reps = []
    for t in range(0, N, step):
        if t in seen:
            continue
        orb = {t}
        stack = [t]
        while stack:
            v = stack.pop()
            for w in [(v * u) % N for u in gens] + [(v + SM) % N, (v - SM) % N]:
                if w not in orb:
                    orb.add(w)
                    stack.append(w)
        seen |= orb
        reps.append(t)
    return reps


def group_tables(N, M, r, t):
    """rp[y] = r^y, S[y] = 1 + r + ... + r^(y-1) (mod N); mul on codes a + N b."""
    rp = [1 % N] * M
    for y in range(1, M):
        rp[y] = rp[y - 1] * r % N
    S = [0] * (M + 1)
    for y in range(1, M + 1):
        S[y] = (S[y - 1] + rp[y - 1]) % N
    return rp, S


def mul(N, M, rp, t, g, h):
    a1, b1, a2, b2 = g % N, g // N, h % N, h // N
    b = b1 + b2
    a = a1 + rp[b1] * a2 + (t if b >= M else 0)
    return a % N + N * (b % M)


def fingerprint(N, M, r, t):
    """vectorised: the order of g = (a, b) is o_b * ord(x^c) with o_b = M / gcd(b, M) and g^o_b = x^c (the y-coordinate
    of g^k is k b mod M, so g^o_b lies in <x>); centre = elements commuting with x and y"""
    rp, S = group_tables(N, M, r, t)
    m = N * M
    g = np.arange(m, dtype=np.int64)
    rpa = np.array(rp, dtype=np.int64)

    def mulv(u, v):
        a1, b1, a2, b2 = u % N, u // N, v % N, v // N
        bb = b1 + b2
        return (a1 + rpa[b1] * a2 + np.where(bb >= M, t, 0)) % N + N * (bb % M)

    ob = M // np.gcd(g // N, M)
    c = np.zeros(m, dtype=np.int64)
    p = g.copy()
    for k in range(1, M + 1):
        sel = ob == k
        c[sel] = p[sel] % N
        if k < M:
            p = mulv(p, g)
    orders = ob * (N // np.gcd(c, N))
    x, y = 1 % N, N % m
    centre = int(np.sum((mulv(g, x) == mulv(x, g)) & (mulv(g, y) == mulv(y, g))))
    cc = math.gcd(N, r - 1) if r % N != 1 else N
    derived = N // cc
    vals, cnts = np.unique(orders, return_counts=True)
    hist = dict(zip(vals.tolist(), cnts.tolist()))
    expo = int(vals.max())
    return f"Z{centre}_D{derived}_e{expo}_" + ",".join(f"{o}:{hist[o]}" for o in sorted(hist)), hist.get(2, 0)


def size_A(N, M, r, t, SM, ls):
    units = [i for i in range(1, N) if math.gcd(i, N) == 1]
    g = math.gcd(SM, N)
    tot = 0
    for i in units:
        for l in ls:
            rhs = (t * (i - l)) % N
            if rhs % g == 0:
                tot += g
    return tot


def main():
    lo, hi, D = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    parity = sys.argv[4]
    allmode = "--all" in sys.argv
    B = l1_ball3(D)
    rows = []
    for n in range(lo, hi + 1):
        if (parity == "odd" and n % 2 == 0) or (parity == "even" and n % 2 == 1):
            continue
        for M in range(2, n + 1):
            if n % M:
                continue
            N = n // M
            if N < 2:
                continue
            facN = factor(N)
            units_M = [l for l in range(1, M) if math.gcd(l, M) == 1]
            seen = set()
            for r in range(1, N):
                if math.gcd(r, N) != 1 or r in seen or pow(r, M, N) != 1:
                    continue
                cls = {pow(r, l, N) for l in units_M}
                seen |= cls
                if r != min(cls):
                    continue
                if r == 1 and not allmode:
                    continue
                c = math.gcd(N, r - 1) if r != 1 else N
                Mc = M * c
                if not allmode and any(r % (p**e) == 1 and M % p for p, e in facN.items()):
                    absorb_split = True
                else:
                    absorb_split = False
                o = mult_order(r, N)
                ls = [l for l in units_M if pow(r, l, N) == r % N]
                SM = sum(pow(r, k, N) for k in range(M)) % N
                for t in t_orbit_reps(N, M, r, c, SM, ls):
                    if t == 0 and absorb_split:
                        continue
                    A = size_A(N, M, r, t, SM, ls)
                    fam_b2a = n**3 / (8 * A)
                    fam_haar = n**3 / (12 * A)
                    fp, ninv = fingerprint(N, M, r, t)
                    sh = ["B2a", "haar"] if Mc <= B else []
                    if ninv >= 1:
                        sh += ["S3", "S1a"]
                    if ninv >= 2:
                        sh += ["S2a", "S2b"]
                    if ninv >= 3:
                        sh += ["S1b", "S1c"]
                    shapes = ",".join(sh) if sh else "-"
                    rows.append((n, N, M, r, t, c, Mc, o, A, fam_b2a, fam_haar, ninv, fp, shapes))
    rows.sort()
    print(f"# {parity}-order metacyclic groups G(N, M, r, t) with |G| in [{lo}, {hi}], D = {D}, lift bound M c <= L3(D) = {B}: "
          f"{len(rows)} presentations, {sum(1 for r in rows if 'B2a' in r[-1])} pass the bound, "
          f"{len({r[-2] for r in rows})} distinct fingerprints; sum fam_B2a (passing) = {sum(r[9] for r in rows if 'B2a' in r[-1]):.3e}, "
          f"sum fam_haar = {sum(r[10] for r in rows if 'B2a' in r[-1]):.3e}; with involutions: {sum(1 for r in rows if r[11] >= 1)}")
    print("# n N M r t c Mc ord |A| fam_B2a fam_haar ninv fp shapes")
    for r in rows:
        print(" ".join(str(v) for v in r[:9]) + f" {r[9]:.4g} {r[10]:.4g} {r[11]} {r[12]} {r[13]}")


if __name__ == "__main__":
    main()
