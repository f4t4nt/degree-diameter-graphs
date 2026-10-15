#!/usr/bin/env python3
"""Reduced candidate prefixes for ANY bitlift base spec (here the 61 three-vertex degree-5 shapes of shapes.py) over the
metacyclic group G(N, M, r, t) = <x, y | x^N = 1, y^M = x^t, y x y^-1 = x^r> (t = 0: the split group Z_N x|_r Z_M).
Element x^a y^b <-> a + N b; (a1,b1)(a2,b2) = (a1 + r^b1 a2 + t [b1 + b2 >= M], (b1 + b2) mod M) as in bitlift.c.

This is gen_b2.py's t-aware sequential-stabiliser scheme (lift_two_vertex/) with the spec passed in instead of looked up and
without the two-vertex swap rule (shape automorphisms are NOT used: sound, over-covering by at most |Aut(shape)| <= 6):
  A = {x -> x^i, y -> x^j y^l : i unit mod N, l unit mod M with r^l = r mod N, j S_M = t (i - l) mod N}
      (S_b = 1 + r + ... + r^(b-1)); phi(x^a y^b) = x^(i a + j S_b + t floor(l b / M)) y^(l b mod M): automorphisms of G
      (checked by brute force in check_k3.py) acting slot-wise on lifts.  A contains Inn(G) = the residual gauge group
      after the spanning tree is fixed to the identity.
  Gamma = A x (inversion of each loop slot) x (permutations of parallel slots: same kind at the same place).
  slot 1 over Gamma-orbit representatives, slot d over orbit representatives of the stabiliser of the earlier slots, a
  parallel slot restricted to type >= type(its predecessor) (strictly larger in (type, index) when the stabiliser is
  trivial), the last free slot emitted as a range [lo, hi) (loops: the lower half of the y-exponents as the inversion
  proxy) or, when a stabiliser survives, as explicit orbit representatives.
Types = Gamma-orbits of G, numbered by their minimum element.  Semi slots take involutions only, loops non-involutions.
Validated by check_k3.py (brute-force orbit coverage on small split and non-split groups).

usage: gen_k3.py N M r t "SPEC" [--count-only] [--out FILE]
output: lines "v1 ... lo hi" (or full lines) on stdout / FILE; "# ... lines, candidates" on stderr.
"""
import sys

import numpy as np


class Group:
    def __init__(self, N, M, r, t):
        self.N, self.M, self.r, self.t = N, M, r, t
        m = self.m = N * M
        rp = [1 % N] * M
        for y in range(1, M):
            rp[y] = rp[y - 1] * r % N
        if rp[M - 1] * r % N != 1 % N:
            raise SystemExit("r^M != 1 mod N")
        if (t * (r - 1)) % N:
            raise SystemExit("t (r - 1) != 0 mod N")
        self.rp = np.array(rp, dtype=np.int64)
        S = [0] * (M + 1)
        for y in range(1, M + 1):
            S[y] = (S[y - 1] + rp[y - 1]) % N
        self.S = np.array(S, dtype=np.int64)
        sM = S[M]
        g = np.arange(m)
        self.X = g % N
        self.Y = g // N
        Yi = (M - self.Y) % M
        self.inv = ((-self.rp[Yi] * (self.X + t * (self.Y != 0))) % N) + N * Yi
        units = [i for i in range(1, N) if np.gcd(i, N) == 1]
        ls = [l for l in range(1, M) if np.gcd(l, M) == 1 and pow(r, l, N) == r % N]
        trip = []
        for i in units:
            for l in ls:
                rhs = (t * (i - l)) % N
                for j in range(N):
                    if (j * sM) % N == rhs:
                        trip.append((i, j, l))
        self.I0 = np.array([q[0] for q in trip], dtype=np.int64)
        self.J0 = np.array([q[1] for q in trip], dtype=np.int64)
        self.L0 = np.array([q[2] for q in trip], dtype=np.int64)
        typ = np.full(m, -1, dtype=np.int64)
        nt = 0
        for h in range(m):
            if typ[h] >= 0:
                continue
            orb = np.concatenate([self.images(self.I0, self.J0, self.L0, h), self.images(self.I0, self.J0, self.L0, int(self.inv[h]))])
            typ[orb] = nt
            nt += 1
        self.typ = typ
        self.ntypes = nt
        self.typ_min = np.array([int(np.min(np.nonzero(typ == q)[0])) for q in range(nt)])
        self.is_invol = (self.inv == g) & (g != 0)
        self.ninv = int(self.is_invol.sum())

    def images(self, I, J, L, h):
        a, b = int(self.X[h]), int(self.Y[h])
        lb = L * b
        return (I * a + J * self.S[b] + self.t * (lb // self.M)) % self.N + self.N * (lb % self.M)

    def mul(self, g, h):
        a1, b1, a2, b2 = g % self.N, g // self.N, h % self.N, h // self.N
        b = b1 + b2
        return (a1 + int(self.rp[b1]) * a2 + (self.t if b >= self.M else 0)) % self.N + self.N * (b % self.M)


def parse_spec(spec):
    slots = []
    for tok in spec.split(";"):
        if tok.startswith("k=") or tok.startswith("deg="):
            continue
        typ, rest = tok.split(":")
        fixed = typ.endswith("*")
        typ = typ.rstrip("*")
        if typ == "link":
            u, v = (int(q) for q in rest.split(","))
        else:
            u = v = int(rest)
        slots.append((typ, u, v, fixed))
    return slots


def generate(G, spec, count_only=False, out=sys.stdout, part=None, split_depth=None, compact=True):
    """returns (lines, candidates); writes the prefix lines unless count_only.
    part = (i, P): only the subtrees whose index at recursion depth split_depth (default min(4, nfree - 2) free slots
    assigned) is i mod P (a deterministic partition of the family into P parts, enumerated in recursion order).
    compact: a last slot with a surviving stabiliser is written as one line "v1 .. v_{nfree-1} -1 h1 h2 ..." (explicit
    orbit representatives) instead of one full line per representative (bitlift.c's compact prefix lines)."""
    m, g = G.m, np.arange(G.m)
    slots = parse_spec(spec)
    free = [s for s, sl in enumerate(slots) if not sl[3]]
    nfree = len(free)
    depth = nfree - 1
    if split_depth is None:
        split_depth = max(1, min(4, nfree - 2))
    leaf = [0]

    def allowed_mask(s):
        kind = slots[s][0]
        if kind == "loop":
            return (g != 0) & (G.inv != g)
        if kind == "semi":
            return G.is_invol.copy()
        if any(q[3] and q[:3] == slots[s][:3] for q in slots):
            return g != 0
        return np.ones(m, dtype=bool)

    pred = {}
    for a, s in enumerate(free):
        for q in free[:a][::-1]:
            if slots[q][:3] == slots[s][:3]:
                pred[s] = q
                break

    def orbit_reps(I, J, L, mask, loop):
        """one representative (the orbit minimum) per orbit of the masked elements under the group (I, J, L) x inversion"""
        if len(I) == 1:
            sel = mask & (g <= G.inv) if loop else mask
            return [int(x) for x in np.nonzero(sel)[0]]
        if len(I) <= 24:
            sel = np.arange(len(I))
        else:
            # a pseudo-random subset generates a subgroup whose orbits refine the true orbits: representatives stay
            # complete (redundancy only)
            sel = np.random.RandomState(len(I) * 7919 + m).choice(len(I), 24, replace=False)
        lb_ = L[sel, None] * G.Y[None, :]
        img = (I[sel, None] * G.X[None, :] + J[sel, None] * G.S[G.Y][None, :] + G.t * (lb_ // G.M)) % G.N + G.N * (lb_ % G.M)
        if loop:
            img = np.vstack([img, G.inv[None, :]])
        rep = np.arange(m)
        while True:
            new = rep.copy()
            for k in range(img.shape[0]):
                new = np.minimum(new, new[img[k]])
                np.minimum.at(new, img[k], new)
            if np.array_equal(new, rep):
                break
            rep = new
        return np.unique(rep[mask]).tolist()

    def stabiliser(I, J, L, h, loop):
        img = G.images(I, J, L, h)
        keep = img == h
        if loop:
            keep |= img == G.inv[h]
        return I[keep], J[keep], L[keep]

    count = [0, 0]
    buf = []

    def flush():
        out.write("".join(buf))
        buf.clear()

    def rec(d, I, J, L, vals):
        if d == split_depth and part is not None:
            idx = leaf[0]
            leaf[0] += 1
            if idx % part[1] != part[0]:
                return
        s = free[d]
        loop = slots[s][0] == "loop"
        trivial = len(I) == 1
        lb = 0
        lo = 0
        extra = None
        if s in pred:
            p = vals[pred[s]]
            tp = int(G.typ[p])
            lb = max(lb, tp)
            if trivial:
                extra = (G.typ > tp) | (g > p)
                lo = min(p + 1, int(G.typ_min[tp + 1]) if tp + 1 < G.ntypes else m)
        if d == depth:
            pre = " ".join(str(vals[q]) for q in free[:depth])
            if not trivial:
                mask = allowed_mask(s)
                if lb > 0:
                    mask &= G.typ >= lb
                reps = orbit_reps(I, J, L, mask, loop)
                if not reps:
                    return
                if not count_only:
                    if compact:
                        buf.append(f"{pre} -1 " + " ".join(str(h) for h in reps) + "\n")
                    else:
                        buf.extend(f"{pre} {h}\n" for h in reps)
                    if len(buf) >= 10000:
                        flush()
                count[0] += 1 if compact else len(reps)
                count[1] += len(reps)
                return
            hi = m
            if loop:
                hi = G.N * (G.M // 2 + 1) if G.M % 2 == 0 else G.N * ((G.M - 1) // 2 + 1)
                lo = max(lo, 1)
            lo = max(lo, int(G.typ_min[lb]) if lb > 0 else 0)
            if hi > lo:
                if not count_only:
                    buf.append(f"{pre} {lo} {hi}\n")
                    if len(buf) >= 10000:
                        flush()
                count[0] += 1
                count[1] += hi - lo
            return
        mask = allowed_mask(s)
        if lb > 0:
            mask &= G.typ >= lb
        if extra is not None:
            mask &= extra
        for h in orbit_reps(I, J, L, mask, loop):
            vals[s] = h
            I2, J2, L2 = stabiliser(I, J, L, h, loop) if not trivial else (I, J, L)
            rec(d + 1, I2, J2, L2, vals)
        vals.pop(s, None)

    rec(0, G.I0, G.J0, G.L0, {})
    flush()
    return count[0], count[1]


def estimate(G, spec):
    """crude reduced-family size: product of the slot ranges over |A|, parallel slots divided by k!, loops by 2"""
    import math
    slots = parse_spec(spec)
    est = 1.0
    groups = {}
    for kind, u, v, fixed in slots:
        if fixed:
            continue
        groups[(kind, u, v)] = groups.get((kind, u, v), 0) + 1
        est *= {"loop": (G.m - 1 - G.ninv) / 2, "semi": G.ninv, "link": G.m}[kind]
    for k in groups.values():
        est /= math.factorial(k)
    return est / len(G.I0)


def main():
    N, M, r, t = (int(v) for v in sys.argv[1:5])
    spec = sys.argv[5]
    count_only = "--count-only" in sys.argv
    outfn = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else None
    part = None
    if "--part" in sys.argv:
        i = sys.argv.index("--part")
        part = (int(sys.argv[i + 1]), int(sys.argv[i + 2]))
    split_depth = int(sys.argv[sys.argv.index("--split-depth") + 1]) if "--split-depth" in sys.argv else None
    compact = "--no-compact" not in sys.argv
    G = Group(N, M, r, t)
    if outfn:
        with open(outfn + ".tmp", "w") as f:
            lines, cands = generate(G, spec, count_only, f, part, split_depth, compact)
        import os
        os.rename(outfn + ".tmp", outfn)
    else:
        lines, cands = generate(G, spec, count_only, sys.stdout, part, split_depth, compact)
    ptxt = f" part={part[0]}/{part[1]}" if part else ""
    print(f"# gen_k3 N={N} M={M} r={r} t={t} spec={spec} |A|={len(G.I0)} ninv={G.ninv} types={G.ntypes} est={estimate(G, spec):.4g}{ptxt}: "
          f"{lines} lines, {cands} candidates (range sizes, before the evaluator's allowed()/simple filter)", file=sys.stderr)


if __name__ == "__main__":
    main()
