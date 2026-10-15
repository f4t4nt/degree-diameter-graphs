"""Fully reduced shape-103 candidate lists over G = Z_N x|_r Z_M (any N, M with r^M = 1 mod N; a kernel is allowed).

Shape 103 = base "k=2;deg=5;loop:0;link*:0,1;link:0,1;link:0,1;loop:1" (free slots a, c2, c3, b in that order), the
shape of the known (5,7) graphs on 5894 and 5896 vertices and of the (5,8) graphs found here.  Element (x, y) <-> x + N*y and
(x1,y1)(x2,y2) = (x1 + r^y1 x2, y1 + y2), as in bitlift --meta N M r; the dart u -> v with voltage g joins (u, h) to
(v, h g).

Symmetry group Gamma of the family (each generator maps a lift to an isomorphic lift; brute-force orbit check in check_b2.py's style):
  A = {x -> x^i, y -> x^j y^t : i a unit mod N, j s_M = 0 mod N, r^t = r mod N, gcd(t, M) = 1}, s_y = 1 + r + ... +
      r^(y-1), acting diagonally on all voltages (automorphisms of G; all of Aut(G) for prime N and faithful action;
      Inn(G) <= A, so the A-orbits are closed under conjugation);
  inversion of a and inversion of b;
  the vertex swap (a, c2, c3, b) -> (b, a, c2^-1, c3^-1);
  the S3 of link permutations realised by switching at vertex 1: (a, c3, c2, b), (a, c2^-1, c3 c2^-1, c2 b c2^-1),
      (a, c3^-1, c2 c3^-1, c3 b c3^-1) and compositions.
Normal form: type(g) = index of the (A x inversion)-orbit of g, numbered by least element.  a is the least element of
its orbit (over the allowed loop voltages: g != e, g != g^-1), ordered by type; (c2, c3) is the lexicographically least
element of its orbit under Stab(a) x S3 (Stab(a) = {phi in A : phi(a) in {a, a^-1}}, orbits marked by brute force);
b runs over inverse-pair representatives (b <= b^-1 as indices) with type(b) >= type(a) (vertex swap; the type is
invariant under A, inversion and conjugation).  Every admissible tuple is Gamma-equivalent to a tuple of this form:
swap the vertices if type(b) < type(a), move a to its orbit representative by A x inv_a, move (c2, c3) to its
Stab(a) x S3 minimum (this conjugates b and may invert a, which inv_a undoes), finally replace b by b^-1 if needed.
Output: lines "a c2 c3 lo hi" (b runs over the indices [lo, hi), maximal runs of the b-mask; the evaluator's
allowed() filter is a no-op on them), count on stderr.
usage: gen103.py N M r [--first lo hi] [--count-only] [--no-degenerate]
  --first lo hi   only a-representatives with index in [lo, hi) (class blocks)
  --no-degenerate skip the (c2, c3) pairs with a, c2, c3 all in <x> when a is in <x> (connectivity then forces b
                  outside <x>; those lines are otherwise written last)
"""
import math
import sys

import numpy as np


def runs(mask):
    """maximal runs of True in a boolean array as (lo, hi) index pairs"""
    d = np.diff(np.concatenate(([0], mask.astype(np.int8), [0])))
    los = np.flatnonzero(d == 1)
    his = np.flatnonzero(d == -1)
    return list(zip(los.tolist(), his.tolist()))


def main():
    N, M, r = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    first = None
    count_only = False
    no_degenerate = False
    args = sys.argv[4:]
    i = 0
    while i < len(args):
        if args[i] == "--first":
            first = (int(args[i + 1]), int(args[i + 2]))
            i += 3
        elif args[i] == "--count-only":
            count_only = True
            i += 1
        elif args[i] == "--no-degenerate":
            no_degenerate = True
            i += 1
        else:
            raise SystemExit(f"unknown arg {args[i]}")
    m = N * M
    rp = [pow(r, y, N) for y in range(M)]
    if rp[M - 1] * r % N != 1 % N or math.gcd(r, N) != 1:
        raise SystemExit("r^M != 1 mod N or r not a unit")
    S = [0] * M
    for y in range(1, M):
        S[y] = (S[y - 1] + rp[y - 1]) % N
    sM = (S[M - 1] + rp[M - 1]) % N
    Sy = np.array(S, dtype=np.int64)
    rpa = np.array(rp, dtype=np.int64)
    g = np.arange(m, dtype=np.int64)
    X, Y = g % N, g // N
    rinv = np.array([rp[(M - y) % M] for y in range(M)], dtype=np.int64)
    INV = ((-rinv[Y] * X) % N) + N * ((M - Y) % M)
    units = [i for i in range(1, N) if math.gcd(i, N) == 1] or [0]
    js = [j for j in range(N) if (j * sM) % N == 0]
    ts = [t for t in range(1, M) if math.gcd(t, M) == 1 and pow(r, t, N) == r % N]
    trip = [(i, j, t) for i in units for j in js for t in ts]
    I0 = np.array([q[0] for q in trip], dtype=np.int64)
    J0 = np.array([q[1] for q in trip], dtype=np.int64)
    T0 = np.array([q[2] for q in trip], dtype=np.int64)

    def images(I, J, T, h):
        return (I * X[h] + J * Sy[Y[h]]) % N + N * ((T * Y[h]) % M)

    def mul(a, b):
        return int((X[a] + rpa[Y[a]] * X[b]) % N + N * ((Y[a] + Y[b]) % M))

    typ = np.full(m, -1, dtype=np.int64)
    nt = 0
    for h in range(m):
        if typ[h] >= 0:
            continue
        orb = np.concatenate([images(I0, J0, T0, h), images(I0, J0, T0, int(INV[h]))])
        typ[orb] = nt
        nt += 1
    loop_ok = (g != 0) & (INV != g)
    is_rep = g <= INV
    inX = Y == 0
    areps = []
    for t in range(nt):
        members = np.flatnonzero(typ == t)
        if loop_ok[members[0]]:
            areps.append(int(members[0]))
    if first is not None:
        areps = [a for a in areps if first[0] <= a < first[1]]
    print(f"# gen103 N={N} M={M} r={r} |G|={m} |A|={len(I0)} (units {len(units)}, j {len(js)}, t {ts}) types={nt} "
          f"a-classes={len(areps)}: {[(a, (int(X[a]), int(Y[a]))) for a in areps]}", file=sys.stderr)

    out = sys.stdout
    buf = []
    nlines = 0
    ncand = 0
    degenerate = []
    marked = np.zeros(m * m, dtype=bool)
    for a in areps:
        img = images(I0, J0, T0, a)
        keep = (img == a) | (img == INV[a])
        Is, Js, Ts = I0[keep], J0[keep], T0[keep]
        bmask = loop_ok & is_rep & (typ >= typ[a])
        branges = runs(bmask)
        bsize = sum(hi - lo for lo, hi in branges)
        marked[:] = False
        marked[:m] = True                      # c2 = e
        marked[::m] = True                     # c3 = e
        marked[g * m + g] = True               # c2 = c3
        nreps = 0
        a_inX = bool(inX[a])
        for c2 in range(1, m):
            row = marked[c2 * m:(c2 + 1) * m]
            for c3 in np.flatnonzero(~row).tolist():
                if row[c3]:
                    continue
                c2i, c3i = int(INV[c2]), int(INV[c3])
                imgs = [(c2, c3), (c3, c2), (c2i, mul(c3, c2i)), (mul(c3, c2i), c2i), (c3i, mul(c2, c3i)), (mul(c2, c3i), c3i)]
                for p2, p3 in imgs:
                    x2 = images(Is, Js, Ts, p2)
                    x3 = images(Is, Js, Ts, p3)
                    marked[x2 * m + x3] = True
                nreps += 1
                if a_inX and inX[c2] and inX[c3]:
                    if not no_degenerate:
                        degenerate.append((a, c2, c3))
                    continue
                for lo, hi in branges:
                    if not count_only:
                        buf.append(f"{a} {c2} {c3} {lo} {hi}\n")
                    nlines += 1
                ncand += bsize
                if len(buf) >= 20000:
                    out.write("".join(buf))
                    buf.clear()
        print(f"#  a={a}={(int(X[a]), int(Y[a]))} type={int(typ[a])} |Stab|={len(Is)} pair-reps={nreps} b-range={bsize} "
              f"({len(branges)} runs) running total {ncand}", file=sys.stderr)
    for a, c2, c3 in degenerate:
        bmask = loop_ok & is_rep & (typ >= typ[a]) & ~inX
        for lo, hi in runs(bmask):
            if not count_only:
                buf.append(f"{a} {c2} {c3} {lo} {hi}\n")
            nlines += 1
            ncand += hi - lo
    out.write("".join(buf))
    print(f"# shape 103 N={N} M={M} r={r}: {ncand} candidates in {nlines} lines ({len(degenerate)} degenerate pair-reps last)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
