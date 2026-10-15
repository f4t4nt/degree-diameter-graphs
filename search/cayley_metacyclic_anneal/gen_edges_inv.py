#!/usr/bin/env python3
"""Stream the edge list of Cay(G, {t} u S u S^-1) on G = Z_N x|_r Z_M to stdout (one "u v" per line, each undirected
edge once), for a connection set with one involution t = x^tx y^(M/2) (a FOUND line of cayley_sa lists the K slot
pairs followed by the involution pair "tx M/2"; pass that list as is).

Group convention: element x^i y^j has index i + N j; (i1, j1)(i2, j2) = (i1 + r^j1 i2 mod N,
j1 + j2 mod M); vertex g is adjacent to g s for s in the connection set.  Edges g ~ g s (s in S) are written once per
generator; the involution's edges g ~ g t are written only when g < g t.

usage: gen_edges_inv.py N M r x1 y1 ... xK yK tx M/2
"""
import sys


def main() -> None:
    nn, mm, r = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rpow = [1] * mm
    for j in range(1, mm):
        rpow[j] = rpow[j - 1] * r % nn
    assert rpow[mm - 1] * r % nn == 1 % nn, "r^M != 1"
    assert mm % 2 == 0, "M must be even"
    vals = list(map(int, sys.argv[4:]))
    assert len(vals) % 2 == 0 and vals[-1] % mm == mm // 2, "last pair must be the involution tx M/2"
    gens = [(vals[2 * t] % nn, vals[2 * t + 1] % mm) for t in range(len(vals) // 2 - 1)]
    tx, half = vals[-2] % nn, mm // 2
    assert tx * (1 + rpow[half]) % nn == 0, "t is not an involution"
    out = sys.stdout
    for j in range(mm):
        base = nn * j
        shifts = [((rpow[j] * sx) % nn, (j + sy) % mm) for sx, sy in gens]
        j2t = (j + half) % mm
        tsh = (rpow[j] * tx) % nn
        lines = []
        for i in range(nn):
            g = base + i
            for sh, j2 in shifts:
                i2 = i + sh
                if i2 >= nn:
                    i2 -= nn
                lines.append(f"{g} {i2 + nn * j2}\n")
            it = i + tsh
            if it >= nn:
                it -= nn
            gt = it + nn * j2t          # g t = x^(i + r^j tx) y^(j + M/2)
            if g < gt:
                lines.append(f"{g} {gt}\n")
        out.write("".join(lines))


if __name__ == "__main__":
    main()
