#!/usr/bin/env python3
"""Stream the edge list of a witness graph on G = Z_N x|_r Z_M to stdout (one "u v" per line).

Group convention: element x^i y^j has index i + N j; (i1, j1)(i2, j2) = (i1 + r^j1 i2 mod N, j1 + j2 mod M).

usage:
  gen_edges.py cayley N M r x1 y1 ... xk yk
      Cay(G, S u S^-1) with S = {x^xi y^yi}: vertex g ~ g s.  Each undirected edge is written once
      (as g -> g s for s in S).
  gen_edges.py lift N M r k c a1x a1y ... akx aky b1x b1y ... bkx bky l1x l1y ... l(c-1)x l(c-1)y
      2-vertex lift: vertex (v, g) has index v * N M + (i + N j).  Edges (0, g) ~ (0, g a_t),
      (1, g) ~ (1, g b_t), (0, g) ~ (1, g l_s) with l_0 = identity.  Each undirected edge written once.
The list is streamed, never held in memory; for n ~ 8e7 and degree 9 it is ~3.6e8 lines (several GB).
"""
import sys


def main() -> None:
    kind = sys.argv[1]
    nn, mm, r = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    rpow = [1] * mm
    for j in range(1, mm):
        rpow[j] = rpow[j - 1] * r % nn
    assert rpow[mm - 1] * r % nn == 1 % nn, "r^M != 1"
    order = nn * mm
    out = sys.stdout
    if kind == "cayley":
        vals = list(map(int, sys.argv[5:]))
        gens = [(vals[2 * t] % nn, vals[2 * t + 1] % mm) for t in range(len(vals) // 2)]
        for j in range(mm):
            base = nn * j
            shifts = [((rpow[j] * sx) % nn, (j + sy) % mm) for sx, sy in gens]
            for i in range(nn):
                g = base + i
                lines = []
                for sh, j2 in shifts:
                    i2 = i + sh
                    if i2 >= nn:
                        i2 -= nn
                    lines.append(f"{g} {i2 + nn * j2}\n")
                out.write("".join(lines))
    elif kind == "lift":
        k, c = int(sys.argv[5]), int(sys.argv[6])
        vals = list(map(int, sys.argv[7:]))
        pairs = [(vals[2 * t] % nn, vals[2 * t + 1] % mm) for t in range(len(vals) // 2)]
        assert len(pairs) == 2 * k + c - 1, "wrong number of voltages"
        loops0, loops1, links = pairs[:k], pairs[k:2 * k], [(0, 0)] + pairs[2 * k:]
        for j in range(mm):
            base = nn * j
            sh0 = [((rpow[j] * sx) % nn, (j + sy) % mm) for sx, sy in loops0]
            sh1 = [((rpow[j] * sx) % nn, (j + sy) % mm) for sx, sy in loops1]
            shl = [((rpow[j] * sx) % nn, (j + sy) % mm) for sx, sy in links]
            for i in range(nn):
                g = base + i
                lines = []
                for sh, j2 in sh0:
                    i2 = i + sh
                    if i2 >= nn:
                        i2 -= nn
                    lines.append(f"{g} {i2 + nn * j2}\n")
                for sh, j2 in shl:
                    i2 = i + sh
                    if i2 >= nn:
                        i2 -= nn
                    lines.append(f"{g} {order + i2 + nn * j2}\n")
                for sh, j2 in sh1:
                    i2 = i + sh
                    if i2 >= nn:
                        i2 -= nn
                    lines.append(f"{order + g} {order + i2 + nn * j2}\n")
                out.write("".join(lines))
    else:
        raise SystemExit("kind must be cayley or lift")


if __name__ == "__main__":
    main()
