"""Reduced candidate list for 6-regular Cayley graphs Cay(G, {a, a^-1, b, b^-1, c, c^-1}) on G = Z_N x|_r Z_M.

Element (x, y) <-> index x + N*y, (x1,y1)(x2,y2) = (x1 + r^y1 x2, y1 + y2) (conventions of lift.c).
Symmetries used: the automorphism group A = <x -> i x (i in Z_N^*), y -> x^j y (j in J)> where
J = {j : j*S_M = 0 mod N}, S_y = 1 + r + ... + r^(y-1); inversion of each generator; permutations of the
three generators.  The level l(g) = min(y, M - y) is invariant under A and inversion, so every orbit of
triples has a representative with l(a) <= l(b) <= l(c).  Enumeration:
  a over A x inv orbit representatives (formula: at level y the A-orbits on x are the classes of
     gcd(x mod t_y, t_y), t_y = gcd(N, j0 S_y), j0 = N / gcd(N, S_M));
  b over Stab_A(a) x inv orbit representatives with l(b) >= l(a), b not in {e, a, a^-1};
  c over elements with l(c) >= l(b), canonical under inversion (y <= (M-1)/2, x <= (N-1)/2 if y = 0),
     and over Stab_A(a, b) x inv orbit representatives when Stab_A(a, b) is nontrivial.
Output lines: "E a b L" (lift.c enumerates c of level >= L itself) and "T a b c" (explicit triples).
Completeness is checked by brute force in check_cay.py (small groups).

usage: gen_cay.py N M r > prefixes.txt       (candidate count on stderr)
"""
import sys
from math import gcd

import numpy as np


def main():
    N, M, r = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    if pow(r, M, N) != 1 or gcd(r, N) != 1 or r % N == 1:
        raise SystemExit("need r^M = 1 mod N, r != 1 (the order of r may be a proper divisor of M)")
    m = N * M
    H = (M - 1) // 2
    rp = [pow(r, y, N) for y in range(M)]
    Sy = [0] * (M + 1)
    for y in range(M):
        Sy[y + 1] = (Sy[y] + rp[y]) % N
    SM = Sy[M]
    j0 = N // gcd(N, SM)
    J = np.arange(0, N, j0, dtype=np.int64)
    units = np.array([i for i in range(1, N) if gcd(i, N) == 1], dtype=np.int64)
    X = np.arange(N, dtype=np.int64)
    out = sys.stdout

    def inv(g):
        x, y = g % N, g // N
        yi = (M - y) % M
        return (-rp[yi] * x) % N + N * yi

    def level(g):
        y = g // N
        return min(y, M - y)

    def t_of(y):
        v = (j0 * Sy[y]) % N
        return gcd(N, v) if v else N

    def arep(g):
        x, y = g % N, g // N
        if y > H:
            g2 = inv(g)
            x, y = g2 % N, g2 // N
        t = t_of(y)
        xm = x % t
        d = gcd(xm, t) if xm else t
        return (d if d < t else 0) + N * y

    areps = sorted({arep(g) for g in range(1, m)})

    SyA = np.array(Sy[:M], dtype=np.int64)

    def stabilizer(a, within=None):
        """Elements (i, j) of A fixing a, as a pair of arrays (I, Jv); within = a group containing them."""
        xa, ya = a % N, a // N
        if within is None:
            Is, Js = [], []
            for i in units:
                mask = ((J * Sy[ya]) % N) == ((1 - int(i)) * xa) % N
                js = J[mask]
                Is.append(np.full(len(js), int(i), dtype=np.int64))
                Js.append(js)
            return np.concatenate(Is), np.concatenate(Js)
        I, Jv = within
        mask = (I * xa + Jv * Sy[ya]) % N == xa
        return I[mask], Jv[mask]

    # inversion as a permutation of the (M, N) element matrix: (x, y) -> (-r^(M-y) x, M - y)
    inv_y = np.array([(M - y) % M for y in range(M)])
    inv_x = np.array([[(-rp[(M - y) % M] * x) % N for x in range(N)] for y in range(M)], dtype=np.int64)
    level_offset = (N * np.arange(M, dtype=np.int64))[:, None]
    row_cache = {}

    def canon_all(stab):
        """stab = (I, Jv), a group of automorphisms.  Canonical index under stab x inv of every element, as an
        array of size m: for level y the images of x are (i x + j S_y) mod N; the translations inside the
        stabiliser form t'_y Z_N, so the minimum over them is a reduction mod t'_y."""
        I, Jv = stab
        trans = Jv[I == 1]
        j1 = int(trans[trans > 0].min()) if (trans > 0).any() else N
        v = (j1 * SyA) % N
        tp = np.array([gcd(N, int(vv)) if vv else N for vv in v], dtype=np.int64)[:, None]
        # one representative j per distinct i; levels with identical (translation parts, t'_y) share a row
        Iu, idx = np.unique(I, return_index=True)
        Ju = Jv[idx]
        C = (Ju[:, None] * SyA[None, :]) % N          # (|Iu|, M)
        best = np.empty((M, N), dtype=np.int64)
        ikey = Iu.tobytes()
        for y in range(M):
            key = (ikey, int(tp[y, 0]), C[:, y].tobytes())
            if key not in row_cache:
                row = np.full(N, N, dtype=np.int64)
                for ii, cc in zip(Iu, C[:, y]):
                    np.minimum(row, ((ii * X + cc) % N) % int(tp[y, 0]), out=row)
                row_cache[key] = row
            best[y] = row_cache[key]
        c1 = best + level_offset
        c2 = best[inv_y[:, None], inv_x] + level_offset[inv_y]
        return np.minimum(c1, c2).reshape(m)

    lev = np.array([level(g) for g in range(m)])
    n_lines = 0
    n_cand = 0

    def ecount(L):
        return (N - 1) // 2 + N * H if L == 0 else N * (H - L + 1)

    ident = np.arange(m)
    for a in areps:
        la = level(a)
        stab_a = stabilizer(a)
        ca = canon_all(stab_a)
        excl = {0, a, inv(a)}
        breps = [int(g) for g in np.nonzero((ca == ident) & (lev >= la))[0] if int(g) not in excl]
        for b in breps:
            lb = level(b)
            stab_ab = stabilizer(b, within=stab_a)
            if len(stab_ab[0]) == 1:
                out.write(f"E {a} {b} {lb}\n")
                n_lines += 1
                n_cand += ecount(lb)
            else:
                cab = canon_all(stab_ab)
                excl2 = excl | {b, inv(b)}
                gy = gcd(gcd(a // N, b // N), M)
                # the y-coordinates must generate Z_M, otherwise Cay(G, S) is disconnected
                creps = [int(g) for g in np.nonzero((cab == ident) & (lev >= lb))[0]
                         if int(g) not in excl2 and gcd(gy, int(g) // N) == 1]
                lines = [f"T {a} {b} {c}\n" for c in creps]
                out.write("".join(lines))
                n_lines += len(lines)
                n_cand += len(lines)
    print(f"# N={N} M={M} r={r} m={m} areps={len(areps)} lines={n_lines} candidates~{n_cand}", file=sys.stderr)


if __name__ == "__main__":
    main()
