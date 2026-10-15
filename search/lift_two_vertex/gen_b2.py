"""Reduced candidate prefixes for the degree-4 two-vertex lift shapes over G(N, M, r, t) = <x, y | x^N = 1, y^M = x^t,
y x y^-1 = x^r> (t = 0: the split group Z_N x|_r Z_M).  Element x^a y^b <-> a + N b; (a1,b1)(a2,b2) =
(a1 + r^b1 a2 + t [b1 + b2 >= M], (b1 + b2) mod M) as in bitlift.c.

Shapes (bitlift base specs; slot order = spec order, the last free slot is left to the evaluator as a range [lo, hi)):
  B2a   k=2;deg=4;loop:0;loop:1;link*:0,1;link:0,1      loop a at 0, loop b at 1, links e and c
  haar  k=2;deg=4;link*:0,1;link:0,1;link:0,1;link:0,1  links e, c1, c2, c3 (the lift is the Haar graph of G)

Reduction (the sequential-stabiliser scheme of gen_prefix.py, made t-aware):
  A = {x -> x^i, y -> x^j y^l : i unit mod N, l unit mod M with r^l = r mod N, j S_M = t (i - l) mod N}
      (S_b = 1 + r + ... + r^(b-1)); phi(x^a y^b) = x^(i a + j S_b + t floor(l b / M)) y^(l b mod M).  Every such map
      is an automorphism of G (relations checked in check_b2.py by brute force) and acts on lifts slot-wise.
  Gamma = A x (inversion of each loop slot).  Types = Gamma-orbits of G, numbered by their minimum element.
  B2a:  a over Gamma-orbit representatives; b over orbit representatives of the setwise stabiliser of {a, a^-1}
        (times inversion of b) with type(b) >= type(a) (vertex swap 0 <-> 1 maps (a, b, c) to (b, a, c^-1));
        c over the whole range [0, m) (c = e is excluded by the evaluator).
  haar: c1 over A-orbit representatives; c2 over Stab(c1)-orbit representatives with type(c2) >= type(c1)
        (and c2 > c1 when the residual stabiliser is trivial); c3 over the range [lo, m) with type(c3) >= type(c2)
        (parallel links are interchangeable).  The vertex swap (all links inverted) and the re-gauging link swap are
        not used (redundancy only, never a loss of coverage).
  semi-edge shapes S3, S2a, S2b, S1a, S1b, S1c (groups with involutions): the same scheme, semi slots restricted to
        involutions, parallel semi-edges at one vertex sorted like parallel links; no vertex-swap rule.
Validated by check_b2.py (brute-force orbits on small split and non-split groups).

usage: gen_b2.py N M r t SHAPE[,SHAPE...] [--count-only] [--files PREFIX]
output: lines "v1 ... lo hi" on stdout (one shape), or with --files PREFIX one file PREFIX_<shape>.txt per shape plus
PREFIX_<shape>.count holding the summary line; "# ... lines, candidates" on stderr.
"""
import sys

import numpy as np

SPECS = {
    "B2a": "k=2;deg=4;loop:0;loop:1;link*:0,1;link:0,1",
    "haar": "k=2;deg=4;link*:0,1;link:0,1;link:0,1;link:0,1",
    # semi-edge shapes (a semi-edge carries an involution); the semi slots come first (few orbit representatives)
    "S3": "k=2;deg=4;semi:0;semi:1;link*:0,1;link:0,1;link:0,1",            # 3 links + a semi-edge at each vertex
    "S2a": "k=2;deg=4;semi:1;semi:1;loop:0;link*:0,1;link:0,1",             # 2 links, loop at 0, two semi-edges at 1
    "S2b": "k=2;deg=4;semi:0;semi:0;semi:1;semi:1;link*:0,1;link:0,1",      # 2 links, two semi-edges at each vertex
    "S1a": "k=2;deg=4;semi:0;semi:1;loop:0;loop:1;link*:0,1",               # 1 link, loop + semi-edge at each vertex
    "S1b": "k=2;deg=4;semi:0;semi:1;semi:1;semi:1;loop:0;link*:0,1",        # 1 link, loop + semi at 0, three semis at 1
    "S1c": "k=2;deg=4;semi:0;semi:0;semi:0;semi:1;semi:1;semi:1;link*:0,1", # 1 link, three semi-edges at each vertex
    # degree-5 two-vertex shapes: L links + loops / semi-edges with L + 2 loops + semis = 5 at each vertex
    "103": "k=2;deg=5;loop:0;link*:0,1;link:0,1;link:0,1;loop:1",           # 3 links, a loop at each vertex (the (5,8) graphs' shape and slot order)
    "H5": "k=2;deg=5;link*:0,1;link:0,1;link:0,1;link:0,1;link:0,1",         # 5 links (Haar graph)
    "L1": "k=2;deg=5;loop:0;loop:0;loop:1;loop:1;link*:0,1",                 # 1 link, two loops at each vertex
    "T4": "k=2;deg=5;semi:0;semi:1;link*:0,1;link:0,1;link:0,1;link:0,1",    # 4 links, a semi-edge at each vertex
    "T3a": "k=2;deg=5;semi:1;semi:1;loop:0;link*:0,1;link:0,1;link:0,1",     # 3 links, loop at 0, two semi-edges at 1
    "T3b": "k=2;deg=5;semi:0;semi:0;semi:1;semi:1;link*:0,1;link:0,1;link:0,1",  # 3 links, two semi-edges at each vertex
    "T2a": "k=2;deg=5;semi:0;semi:1;loop:0;loop:1;link*:0,1;link:0,1",       # 2 links, loop + semi-edge at each vertex
    "T2b": "k=2;deg=5;semi:0;semi:1;semi:1;semi:1;loop:0;link*:0,1;link:0,1",  # 2 links, loop + semi at 0, three semis at 1
    "T2c": "k=2;deg=5;semi:0;semi:0;semi:0;semi:1;semi:1;semi:1;link*:0,1;link:0,1",  # 2 links, three semi-edges at each vertex
    "T1b": "k=2;deg=5;semi:1;semi:1;loop:0;loop:0;loop:1;link*:0,1",         # 1 link, two loops at 0, loop + two semis at 1
    "T1c": "k=2;deg=5;semi:1;semi:1;semi:1;semi:1;loop:0;loop:0;link*:0,1",  # 1 link, two loops at 0, four semis at 1
    "T1d": "k=2;deg=5;semi:0;semi:0;semi:1;semi:1;loop:0;loop:1;link*:0,1",  # 1 link, loop + two semis at each vertex
    "T1e": "k=2;deg=5;semi:0;semi:0;semi:1;semi:1;semi:1;semi:1;loop:0;link*:0,1",  # 1 link, loop + two semis at 0, four semis at 1
    "T1f": "k=2;deg=5;semi:0;semi:0;semi:0;semi:0;semi:1;semi:1;semi:1;semi:1;link*:0,1",  # 1 link, four semi-edges at each vertex
}


def symmetric_rule(slots):
    """vertex-swap rule: if the non-link slot kinds at the two vertices agree (as multisets), the lift with the vertices
    swapped (links inverted) is isomorphic, so WLOG type(first free non-link slot at 1) >= type(first at 0); returns
    {slot at vertex 1: slot at vertex 0} or {} (the first slots must be of the same kind: specs list semis before loops)"""
    k0 = sorted(s[0] for s in slots if s[0] != "link" and s[1] == 0)
    k1 = sorted(s[0] for s in slots if s[0] != "link" and s[1] == 1)
    if not k0 or k0 != k1:
        return {}
    f0 = next(i for i, s in enumerate(slots) if s[0] != "link" and s[1] == 0 and not s[3])
    f1 = next(i for i, s in enumerate(slots) if s[0] != "link" and s[1] == 1 and not s[3])
    if slots[f0][0] != slots[f1][0]:
        return {}
    return {f1: f0}


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
        # inverse of (x, y): y' = (M - y) mod M, x' = -r^y' (x + t [y != 0])
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


def generate(G, shape, count_only=False, out=sys.stdout):
    """returns (lines, candidates); writes the prefix lines unless count_only"""
    m, g = G.m, np.arange(G.m)
    slots = parse_spec(SPECS[shape])
    free = [s for s, sl in enumerate(slots) if not sl[3]]
    nfree = len(free)
    depth = nfree - 1

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
    # vertex-swap rule for the symmetric shapes (B2a, 103, L1, S3, ...): type(first slot at 1) >= type(first slot at 0)
    vrule = symmetric_rule(slots)

    def orbit_reps(I, J, L, mask, loop):
        """one representative (the orbit minimum) per orbit of the masked elements under the group (I, J, L) x inversion"""
        if len(I) == 1:
            sel = mask & (g <= G.inv) if loop else mask
            return [int(x) for x in np.nonzero(sel)[0]]
        # orbits = connected components of the maps h -> phi_k(h) for a generating subset of the stabiliser (all its
        # elements when it is small, otherwise 24 pseudo-random ones: they generate a subgroup, whose orbits refine
        # the true orbits, so the representatives stay complete; redundancy only) plus inversion for loops;
        # minimum-label propagation on the m elements gives the orbit minimum of every element
        if len(I) <= 24:
            sel = np.arange(len(I))
        else:
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
        s = free[d]
        loop = slots[s][0] == "loop"
        trivial = len(I) == 1
        lb = 0
        if s in vrule:
            lb = int(G.typ[vals[vrule[s]]])
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
                # the residual stabiliser of the earlier slots still acts on this slot (e.g. loops inside <x> are fixed
                # by every y -> x^j y): emit its orbit representatives as full lines instead of a range
                mask = allowed_mask(s)
                if lb > 0:
                    mask &= G.typ >= lb
                reps = orbit_reps(I, J, L, mask, loop)
                if not count_only:
                    buf.extend(f"{pre} {h}\n" for h in reps)
                    if len(buf) >= 10000:
                        flush()
                count[0] += len(reps)
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


def main():
    N, M, r, t = (int(v) for v in sys.argv[1:5])
    shapes = sys.argv[5].split(",")
    count_only = "--count-only" in sys.argv
    files = sys.argv[sys.argv.index("--files") + 1] if "--files" in sys.argv else None
    G = Group(N, M, r, t)
    for shape in shapes:
        if files:
            with open(f"{files}_{shape}.txt.tmp", "w") as f:
                lines, cands = generate(G, shape, count_only, f)
        else:
            lines, cands = generate(G, shape, count_only)
        summary = (f"# gen_b2 N={N} M={M} r={r} t={t} shape={shape} |A|={len(G.I0)} types={G.ntypes}: {lines} lines, "
                   f"{cands} candidates (range sizes, before the evaluator's allowed()/simple filter)")
        print(summary, file=sys.stderr)
        if files:
            import os
            os.rename(f"{files}_{shape}.txt.tmp", f"{files}_{shape}.txt")
            with open(f"{files}_{shape}.count", "w") as f:
                f.write(summary + "\n")


if __name__ == "__main__":
    main()
