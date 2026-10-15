#!/usr/bin/env python3
"""Independent edge-list builder for a 2-vertex voltage lift over G(N, M, r, t) = <x, y | x^N = 1, y^M = x^t,
y x y^-1 = x^r> (own implementation of the group law; shares nothing with bitlift.c / gen_b2.py except the conventions):
element x^a y^b <-> a + N b, lift vertex (u, g) <-> u*|G| + g, dart u -> v with voltage w joins (u, g) to (v, g w);
a loop at u with voltage a gives (u, g) ~ (u, g a) (both directions arise), a semi-edge with involution s gives
(u, g) ~ (u, g s), a link u -> v with voltage w gives (u, g) ~ (v, g w).
usage: make_edges.py N M r t "SPEC" "[v0,v1,...]" OUT.edges     (SPEC = bitlift base spec, voltages in slot order)
Prints order, edge count and the degree range; exits 1 if the lift is not a simple regular graph of the spec's degree.
"""
import sys


def main():
    N, M, r, t = (int(v) for v in sys.argv[1:5])
    spec, volts, out = sys.argv[5], sys.argv[6], sys.argv[7]
    vals = [int(v) for v in volts.strip("[] ").split(",")]
    assert pow(r, M, N) == 1 % N and (t * (r - 1)) % N == 0, "not a group presentation"
    m = N * M
    rp = [pow(r, b, N) for b in range(M)]

    def mul(g, h):
        a1, b1, a2, b2 = g % N, g // N, h % N, h // N
        b = b1 + b2
        return (a1 + rp[b1] * a2 + (t if b >= M else 0)) % N + N * (b % M)

    slots = []
    deg = None
    for tok in spec.split(";"):
        if tok.startswith("k="):
            assert int(tok[2:]) == 2
            continue
        if tok.startswith("deg="):
            deg = int(tok[4:])
            continue
        kind, rest = tok.split(":")
        kind = kind.rstrip("*")
        if kind == "link":
            u, v = (int(q) for q in rest.split(","))
        else:
            u = v = int(rest)
        slots.append((kind, u, v))
    assert len(slots) == len(vals), "voltage count != slot count"
    edges = set()
    for g in range(m):
        for (kind, u, v), w in zip(slots, vals):
            h = mul(g, w)
            a, b = u * m + g, v * m + h
            if a == b:
                raise SystemExit("loop edge: not simple")
            edges.add((min(a, b), max(a, b)))
    degs = [0] * (2 * m)
    for a, b in edges:
        degs[a] += 1
        degs[b] += 1
    with open(out, "w") as f:
        for a, b in sorted(edges):
            f.write(f"{a} {b}\n")
    print(f"order {2 * m} edges {len(edges)} degree min {min(degs)} max {max(degs)} (spec degree {deg})")
    if min(degs) != deg or max(degs) != deg:
        raise SystemExit("not regular of the spec degree (multi-edge collapsed)")


if __name__ == "__main__":
    main()
