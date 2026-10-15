"""Reduced candidate prefixes for an arbitrary base pregraph over G = Z_N x|_r Z_M (r^M = 1 mod N; a kernel is allowed).

Sequential-stabiliser scheme: Gamma = A x (inversion of every loop slot), where
A = {x -> x^i, y -> x^j y^t : i a unit mod N, j s_M = 0 mod N, r^t = r, gcd(t, M) = 1} (s_M = 1 + r + ... + r^(M-1))
is a group of automorphisms of G (all of Aut(G) when N is prime and the action is faithful; t != 1 only with a kernel).  Element (x, y) <-> x + N*y as in
bitlift --meta N M r.  Free slots are taken in spec order: slot 1 over Gamma-orbit representatives of its allowed
values, slot 2 over orbit representatives of the setwise stabiliser of slot 1's value ({a, a^-1} for a loop), and so
on; the last free slot is left to the evaluator as an index range [lo, hi).

Base-graph symmetries (all type-preserving, type = Gamma-orbit index of the element, which conjugation, A and
inversion preserve):
  * parallel slots of the same kind at the same place (loops at one vertex, semi-edges at one vertex, non-tree links
    between the same pair): type-sorted; index-sorted (strict) when the residual stabiliser is trivial;
  * vertex chain rules (--rules JSON "[[u0, [others]], ...]"): T(u0) <=lex T(u) for the listed u,
    T(u) = sorted tuple of loop types at vertex u (sound: the vertex automorphism group acts transitively on each
    rule set inside the stabiliser of the previous chain vertices).
Validated by brute-force orbit enumeration on small groups (check_b2.py does the same for the shapes of gen_b2.py).

usage: gen_prefix.py N M r "SPEC" [--rules JSON] [--depth d] [--count-only] [--first lo hi]
                     [--loop-minorder o] [--loop-exclude-levels j1,j2,...]
output: lines "v_1 ... v_{d} lo hi" (stdout; d = nfree - 1 by default), candidate count on stderr.

Heuristic strata restriction (not a complete sweep when used): --loop-minorder o excludes loop voltages
of order < o (a loop of order o puts o-cycles into the lift), --loop-exclude-levels excludes loop voltages x^i y^j with
min(j, M - j) in the list (level 0 = <x>, level M/2 = z <x> are the hopeless layers of the B2a sweeps).  Only loop
slots that are generated here are filtered; a loop left to the evaluator as the last slot is not.
"""
import json
import sys

import numpy as np


def parse_spec(spec):
    slots = []
    k = deg = None
    for tok in spec.split(";"):
        if tok.startswith("k="):
            k = int(tok[2:])
            continue
        if tok.startswith("deg="):
            deg = int(tok[4:])
            continue
        typ, rest = tok.split(":")
        fixed = typ.endswith("*")
        typ = typ.rstrip("*")
        if typ == "link":
            u, v = (int(t) for t in rest.split(","))
        else:
            u = v = int(rest)
        slots.append((typ, u, v, fixed))
    return k, deg, slots


class Group:
    """Z_N x|_r Z_M with the automorphism family A and element types."""

    def __init__(self, N, M, r):
        self.N, self.M, self.r = N, M, r
        m = self.m = N * M
        rp = [1] * M
        for y in range(1, M):
            rp[y] = rp[y - 1] * r % N
        if rp[M - 1] * r % N != 1 % N:
            raise SystemExit(f"r^M != 1 mod N")
        self.rp = rp
        S = [0] * M  # S[y] = 1 + r + ... + r^(y-1)
        for y in range(1, M):
            S[y] = (S[y - 1] + rp[y - 1]) % N
        sM = (S[M - 1] + rp[M - 1]) % N
        self.S = np.array(S, dtype=np.int64)
        g = np.arange(m)
        self.X = g % N
        self.Y = g // N
        rinv = np.array([rp[(M - y) % M] for y in range(M)], dtype=np.int64)
        self.inv = ((-rinv[self.Y] * self.X) % N) + N * ((M - self.Y) % M)
        self.e = 0
        units = [i for i in range(1, N) if np.gcd(i, N) == 1]
        js = [j for j in range(N) if (j * sM) % N == 0]
        # y -> y^t is an automorphism when r^t = r (same action) and gcd(t, M) = 1 (non-trivial only with a kernel)
        ts = [t for t in range(1, M) if np.gcd(t, M) == 1 and pow(r, t, N) == r % N]
        trip = [(i, j, t) for i in units for j in js for t in ts]
        self.I0 = np.array([q[0] for q in trip], dtype=np.int64)
        self.J0 = np.array([q[1] for q in trip], dtype=np.int64)
        self.T0 = np.array([q[2] for q in trip], dtype=np.int64)
        # element types: orbits of A x inversion, numbered by minimum element
        typ = np.full(m, -1, dtype=np.int64)
        nt = 0
        for h in range(m):
            if typ[h] >= 0:
                continue
            orb = np.concatenate([self.images(self.I0, self.J0, self.T0, h), self.images(self.I0, self.J0, self.T0, int(self.inv[h]))])
            typ[orb] = nt
            nt += 1
        self.typ = typ
        self.ntypes = nt
        self.typ_min = np.array([int(np.min(np.nonzero(typ == t)[0])) for t in range(nt)])
        self.is_invol = (self.inv == g) & (g != self.e)
        # element orders (vectorised repeated multiplication; the maximum order is at most |G|)
        order = np.zeros(m, dtype=np.int64)
        rpa = np.array(rp, dtype=np.int64)
        px, py = self.X.copy(), self.Y.copy()      # p = g^t, as (x, y)
        t = 1
        while True:
            done = (px == 0) & (py == 0) & (order == 0)
            order[done] = t
            if np.all(order > 0):
                break
            px, py = (px + rpa[py] * self.X) % N, (py + self.Y) % M   # p * g
            t += 1
            if t > m:
                raise SystemExit("order computation failed")
        self.order = order
        self.level = np.minimum(self.Y, M - self.Y)

    def images(self, I, J, T, h):
        return (I * self.X[h] + J * self.S[self.Y[h]]) % self.N + self.N * ((T * self.Y[h]) % self.M)

    def mul(self, g, h):
        return ((self.X[g] + self.rp[self.Y[g]] * self.X[h]) % self.N) + self.N * ((self.Y[g] + self.Y[h]) % self.M)


def main():
    N, M, r = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    spec = sys.argv[4]
    rules = []
    depth = None
    count_only = False
    first_range = None  # --first lo hi: only first-slot representatives with index in [lo, hi) (class blocks)
    loop_minorder = 0
    loop_excl_levels = []
    slot_orders = {}
    args = sys.argv[5:]
    i = 0
    while i < len(args):
        if args[i] == "--rules":
            rules = json.loads(args[i + 1])
            i += 2
        elif args[i] == "--depth":
            depth = int(args[i + 1])
            i += 2
        elif args[i] == "--first":
            first_range = (int(args[i + 1]), int(args[i + 2]))
            i += 3
        elif args[i] == "--loop-minorder":
            loop_minorder = int(args[i + 1])
            i += 2
        elif args[i] == "--loop-exclude-levels":
            loop_excl_levels = [int(t) for t in args[i + 1].split(",") if t.strip()]
            i += 2
        elif args[i] == "--slot-orders":
            # SLOT:LO:HI (repeatable): the loop slot SLOT (spec index) only takes voltages of order in [LO, HI]
            s_, lo_, hi_ = (int(t) for t in args[i + 1].split(":"))
            slot_orders[s_] = (lo_, hi_)
            i += 2
        elif args[i] == "--count-only":
            count_only = True
            i += 1
        else:
            raise SystemExit(f"unknown arg {args[i]}")
    G = Group(N, M, r)
    m, g = G.m, np.arange(G.m)
    k, deg, slots = parse_spec(spec)
    free = [s for s, sl in enumerate(slots) if not sl[3]]
    nfree = len(free)
    if depth is None:
        depth = nfree - 1
    if depth < 1 or depth > nfree:
        raise SystemExit("bad depth")
    loop_ok = np.ones(m, dtype=bool)
    if loop_minorder:
        loop_ok &= G.order >= loop_minorder
    for lv in loop_excl_levels:
        loop_ok &= G.level != lv
    if loop_minorder or loop_excl_levels:
        adm = (g != G.e) & (G.inv != g)
        print(f"# loop filter: minorder {loop_minorder}, excluded levels {loop_excl_levels}: "
              f"{int((adm & loop_ok).sum())} of {int(adm.sum())} loop-admissible elements kept", file=sys.stderr)

    def allowed_mask(s):
        typ = slots[s][0]
        if typ == "loop":
            ok = (g != G.e) & (G.inv != g) & loop_ok
            if s in slot_orders:
                ok &= (G.order >= slot_orders[s][0]) & (G.order <= slot_orders[s][1])
            return ok
        if typ == "semi":
            return G.is_invol.copy()
        # a non-tree link may carry the identity unless a tree link joins the same pair (double edge)
        if any(t[3] and t[:3] == slots[s][:3] for t in slots):
            return g != G.e
        return np.ones(m, dtype=bool)

    # predecessor of the same kind at the same place (for sorting), per free slot
    pred = {}
    for a, s in enumerate(free):
        for t in free[:a][::-1]:
            if slots[t][:3] == slots[s][:3]:
                pred[s] = t
                break
    # loop slots per vertex (in spec order) and the vertex chain rules
    loops_at = {}
    for s in free:
        if slots[s][0] == "loop":
            loops_at.setdefault(slots[s][1], []).append(s)
    refs = {}  # vertex u -> list of reference vertices u0 with T(u0) <= T(u)
    for u0, others in rules:
        for u in others:
            refs.setdefault(u, []).append(u0)
    for u, lst in refs.items():
        for u0 in lst:
            if u0 not in loops_at or u not in loops_at or len(loops_at[u0]) != len(loops_at[u]):
                raise SystemExit("rule vertices must carry the same number of loops")
            if max(loops_at[u0]) > min(loops_at[u]):
                raise SystemExit("rule reference vertex must precede in slot order")

    def orbit_reps(I, J, T, mask, loop):
        if len(I) == 1:
            sel = mask & (g <= G.inv) if loop else mask
            return [int(x) for x in np.nonzero(sel)[0]]
        reps = []
        seen = ~mask
        for h in range(m):
            if seen[h]:
                continue
            seen[G.images(I, J, T, h)] = True
            if loop:
                seen[G.images(I, J, T, int(G.inv[h]))] = True
            reps.append(h)
        return reps

    def stabiliser(I, J, T, h, loop):
        img = G.images(I, J, T, h)
        keep = img == h
        if loop:
            keep |= img == G.inv[h]
        return I[keep], J[keep], T[keep]

    def lower_type_bound(s, vals):
        """minimum type allowed for loop slot s by the vertex chain rules, given the values so far"""
        u = slots[s][1]
        lb = 0
        if slots[s][0] == "loop" and u in refs:
            j = loops_at[u].index(s)
            partial = [int(G.typ[vals[t]]) for t in loops_at[u][:j]]
            for u0 in refs[u]:
                T0 = sorted(int(G.typ[vals[t]]) for t in loops_at[u0])
                if partial == T0[:j]:
                    lb = max(lb, T0[j])
                elif partial < T0[:j]:
                    return None  # cannot happen (the partial tuple is never allowed below T0)
        return lb

    out = sys.stdout
    count = [0, 0]  # lines, candidates (sum of range sizes)
    buf = []

    def flush():
        out.write("".join(buf))
        buf.clear()

    def rec(d, I, J, T, vals):
        s = free[d]
        kind = slots[s][0]
        loop = kind == "loop"
        trivial = len(I) == 1
        lb = lower_type_bound(s, vals)
        if lb is None:
            return
        # parallel slot of the same kind at the same place: sort by (type, index); the index part only with a
        # trivial residual stabiliser (the orbit representatives are then exactly the inverse-pair representatives)
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
            # the remaining slots (normally just the last one) go to the evaluator; emit the range for slot s
            hi = m
            if loop:
                hi = N * (M // 2 + 1) if M % 2 == 0 else N * ((M - 1) // 2 + 1)
                lo = max(lo, 1)
            lo = max(lo, int(G.typ_min[lb]) if lb > 0 else 0)
            if hi > lo:
                pre = " ".join(str(vals[t]) for t in free[:depth])
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
        reps = orbit_reps(I, J, T, mask, loop)
        if d == 0 and first_range is not None:
            reps = [h for h in reps if first_range[0] <= h < first_range[1]]
        for h in reps:
            vals[s] = h
            if d + 1 == depth and depth == nfree:
                pre = " ".join(str(vals[t]) for t in free)
                if not count_only:
                    buf.append(pre + "\n")
                    if len(buf) >= 10000:
                        flush()
                count[0] += 1
                count[1] += 1
            else:
                I2, J2, T2 = stabiliser(I, J, T, h, loop) if not trivial else (I, J, T)
                rec(d + 1, I2, J2, T2, vals)
        vals.pop(s, None)

    rec(0, G.I0, G.J0, G.T0, {})
    flush()
    print(f"# gen_prefix N={N} M={M} r={r} |A|={len(G.I0)} types={G.ntypes} free={free} depth={depth}: "
          f"{count[0]} lines, {count[1]} candidates (range sizes, before the evaluator's allowed() filter)", file=sys.stderr)


if __name__ == "__main__":
    main()
