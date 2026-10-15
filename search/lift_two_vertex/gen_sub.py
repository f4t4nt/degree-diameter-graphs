"""Sub-block (i, j) of the S_3-reduced shape-103 family over G = Z_N x|_r Z_M in gen103.py's normal form, with O(m)
memory (gen103.py marks orbits in an m x m byte array, 0.9-1.4 GB at |G| = 3e4; this script canonicalises instead).

Shape 103 = "k=2;deg=5;loop:0;link*:0,1;link:0,1;link:0,1;loop:1", element (x, y) <-> x + N y, group law
(x1, y1)(x2, y2) = (x1 + r^y1 x2, y1 + y2).  Sub-block (i, j), 1 <= i <= j <= (M-1)/2, gcd(i, M) = gcd(j, M) = 1:
loop a = y^i at base vertex 0 (the least element of its A x inversion class), loop b = x^s y^j at base vertex 1 for
every s (layer j is a single A-class with type >= type(a), all its elements are inverse-pair representatives), links
(c2, c3) = the lexicographically least element of the Stab(a) x S_3 orbit, exactly as in gen103.py.  Stab(a) is
computed by brute force over A = {x -> x^u, y -> x^t y^w} and must act on x-coordinates only (x -> u x, u in a
subgroup U of the units); then the least element of a Stab(a)-orbit of a pair is computed in closed form:
  min_u index(u . p2) = N y2 + cmin[x2]   (cmin[x] = min of the U-orbit of x),
  if x2 != 0 the minimising u is unique (free action of units on Z_N^*, N prime), which fixes p3's image;
  if x2 == 0 all u in U tie and p3 is minimised independently.
The S_3 images are gen103's: (c2, c3), (c3, c2), (c2^-1, c3 c2^-1), (c3 c2^-1, c2^-1), (c3^-1, c2 c3^-1), (c2 c3^-1, c3^-1).
Output: lines "a c2 c3 lo hi" with [lo, hi) = [j N, (j+1) N); the count on stderr.  Validate against gen103.py with
  gen103.py N M r --first iN iN+1 | awk '{print $1, $2, $3}' | sort -u   versus   gen_sub.py N M r i j | awk '{print $1, $2, $3}' | sort -u
usage: gen_sub.py N M r i j
"""
import math
import sys

import numpy as np


def main():
    N, M, r, i, j = (int(v) for v in sys.argv[1:6])
    if not (1 <= i <= j <= (M - 1) // 2) or math.gcd(i, M) != 1 or math.gcd(j, M) != 1:
        raise SystemExit("need 1 <= i <= j <= (M-1)/2 with i, j coprime to M")
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
    units = [u for u in range(1, N) if math.gcd(u, N) == 1]
    js = [t for t in range(N) if (t * sM) % N == 0]
    ts = [w for w in range(1, M) if math.gcd(w, M) == 1 and pow(r, w, N) == r % N]
    trip = [(u, t, w) for u in units for t in js for w in ts]
    I0 = np.array([q[0] for q in trip], dtype=np.int64)
    J0 = np.array([q[1] for q in trip], dtype=np.int64)
    T0 = np.array([q[2] for q in trip], dtype=np.int64)

    def images(I, J, T, h):
        return (I * X[h] + J * Sy[Y[h]]) % N + N * ((T * Y[h]) % M)

    a = N * i
    # a must be the least element of its (A x inv)-class (gen103's class representative)
    orb = np.concatenate([images(I0, J0, T0, a), images(I0, J0, T0, int(INV[a]))])
    if int(orb.min()) != a:
        raise SystemExit(f"a = y^{i} is not the least element of its class (min {int(orb.min())})")
    img = images(I0, J0, T0, a)
    keep = (img == a) | (img == INV[a])
    Us, Js, Ts = I0[keep], J0[keep], T0[keep]
    if np.any(Js != 0) or np.any(Ts != 1):
        raise SystemExit("Stab(a) does not act on x-coordinates only; use gen103.py")
    U = sorted(set(int(u) for u in Us))
    if len(U) != len(units):
        print(f"# note: Stab(a) is a proper subgroup of the units (order {len(U)})", file=sys.stderr)
    # cmin[x] = min of the U-orbit of x; umin[x] = the u with u x = cmin[x] (unique for x != 0 when N is prime)
    Ua = np.array(U, dtype=np.int64)
    cmin = np.zeros(N, dtype=np.int64)
    umin = np.ones(N, dtype=np.int64)
    for x in range(1, N):
        vals = (Ua * x) % N
        k = int(np.argmin(vals))
        cmin[x] = vals[k]
        umin[x] = Ua[k]
        if np.sum(vals == vals[k]) != 1:
            raise SystemExit("non-free action of U on Z_N^*: N must be prime")
    print(f"# gen_sub N={N} M={M} r={r} |G|={m} |A|={len(I0)} |Stab(a)|={len(U)} a=y^{i} (index {a}) b-layer {j} "
          f"[{j * N}, {(j + 1) * N})", file=sys.stderr)

    def mul(p, q):
        return (X[p] + rpa[Y[p]] * X[q]) % N + N * ((Y[p] + Y[q]) % M)

    def canon(p2, p3):
        """lexicographically least (index pair, as key p2 * m + p3) of the U-orbit of the pairs (p2, p3)"""
        x2, y2, x3, y3 = X[p2], Y[p2], X[p3], Y[p3]
        q2 = N * y2 + cmin[x2]
        # x2 != 0: unique u; x2 == 0: u free, minimise p3
        q3 = np.where(x2 != 0, N * y3 + (umin[x2] * x3) % N, N * y3 + cmin[x3])
        return q2 * m + q3

    lo, hi = j * N, (j + 1) * N
    out = sys.stdout
    buf = []
    nreps = 0
    c3all = np.arange(m, dtype=np.int64)
    for c2 in range(1, m):
        valid = (c3all != 0) & (c3all != c2)
        c3 = c3all[valid]
        c2a = np.full(c3.shape, c2, dtype=np.int64)
        c2i = np.full(c3.shape, int(INV[c2]), dtype=np.int64)
        c3i = INV[c3]
        m32 = mul(c3, c2i)
        m23 = mul(c2a, c3i)
        key = canon(c2a, c3)
        key = np.minimum(key, canon(c3, c2a))
        key = np.minimum(key, canon(c2i, m32))
        key = np.minimum(key, canon(m32, c2i))
        key = np.minimum(key, canon(c3i, m23))
        key = np.minimum(key, canon(m23, c3i))
        reps = c3[key == c2 * m + c3]
        for c in reps.tolist():
            buf.append(f"{a} {c2} {c} {lo} {hi}\n")
        nreps += len(reps)
        if len(buf) >= 20000:
            out.write("".join(buf))
            buf.clear()
    out.write("".join(buf))
    print(f"# sub-block ({i}, {j}) of N={N} M={M} r={r}: {nreps} pair-reps x {N} b = {nreps * N} candidates", file=sys.stderr)


if __name__ == "__main__":
    main()
