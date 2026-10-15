#!/usr/bin/env python3
"""classcount.py: the number of conjugacy classes k(G) and the commuting probability cp = k(G)/|G| of a metacyclic group
G = Z_N x|_r Z_M = <x, y | x^N = y^M = 1, y x y^-1 = x^r>  (r^M = 1 mod N; faithful or not).

Derivation.  Z_N is normal with abelian quotient Z_M, so conjugation preserves the y-exponent and every class lies in one
layer L_b = {x^a y^b : a in Z_N}.  For g = x^c y^e,
    g (x^a y^b) g^-1 = x^{r^e a + c (1 - r^b)} y^b,
so G acts on L_b = Z_N through the affine maps a -> r^e a + c (1 - r^b).  The translations c (1 - r^b) form the subgroup
d_b Z_N with d_b = gcd(N, r^b - 1); their orbits are the residues a mod d_b, and the multiplications by r (well defined on
Z_{d_b} because r^b = 1 mod d_b) act on them.  Hence
    #classes in layer b = #orbits of <r> on Z_{d_b}
                        = (1/o_b) sum_{k=0}^{o_b - 1} gcd(d_b, r^k - 1)          (Burnside, o_b = ord of r mod d_b)
                        = sum_{e | d_b} phi(e) / ord_e(r)                       (orbits on the elements of exact order e),
and k(G) = sum_{b in Z_M} #classes in layer b.  Exchanging the sums (e | d_b = gcd(N, r^b - 1) iff ord_e(r) | b, which
happens for M / ord_e(r) values of b in Z_M) gives the closed form

    k(G) = M * sum_{e | N} phi(e) / ord_e(r)^2,        cp = k(G) / |G| = sum_{e | N} (phi(e) / N) / ord_e(r)^2,

i.e. cp is the mean of 1 / ord_e(r)^2 over the order e of a uniformly random element of Z_N.  Equivalently, since x^a y^b
and x^c y^d commute iff a (1 - r^d) = c (1 - r^b) mod N, which for random a, c has probability gcd(d_b, d_d) / N,
    cp = (1 / (N M^2)) sum_{b, d in Z_M} gcd(d_b, d_d):
the commuting probability counts the pairs of layers whose elements commute with high probability, i.e. the layer pairs
that produce commutator relators [u, v] = 1 between words of a Cayley graph.  The term e = N is
phi(N) / (N ord_N(r)^2) ~ 1 / M^2 for a faithful presentation, e = 1 gives 1 / N, and every proper divisor e with a small
acting order ord_e(r) adds (phi(e) / N) (M / ord_e(r))^2 / M^2: this is the whole structural content of cp beyond M.
Checks: r = 1 gives N M (abelian); the dihedral Z_N x|_{-1} Z_2 (N odd) gives 1 + (N-1)/2 + 1 = (N+3)/2; N prime and
ord r = M give (N-1)/M + M.  For the one-involution groups (N odd composite, M even, ord r = M, s = r^{M/2}, p = gcd(N, 1+s)
involutions, q = N/p) the divisors of q have ord_e(r) | M/2, so q contributes >= (phi(q)/N) (M/o_q)^2 / M^2 ~ 4/(p M^2)
or more: the acting orders o_p = ord_p(r), o_q = ord_q(r) are what moves cp at fixed (N, M).

`python3 classcount.py N M r` prints k(G), cp and the per-layer counts; `python3 classcount.py test` runs the brute-force
validation (orbit enumeration of the conjugation action on all |G| elements) on small groups incl. even M, composite N,
non-faithful r and r = 1, and the two closed formulas against each other.
The (5,7) group choice used cp as the selection rule: groups with the lowest cp in a window are the softest."""
import math
import sys


def divisors(n):
    small, large = [], []
    i = 1
    while i * i <= n:
        if n % i == 0:
            small.append(i)
            if i * i != n:
                large.append(n // i)
        i += 1
    return small + large[::-1]


def order_mod(r, e, M):
    """multiplicative order of r modulo e, given that it divides M (r^M = 1 mod N, e | N)."""
    if e == 1:
        return 1
    r %= e
    for d in divisors(M):
        if pow(r, d, e) == 1:
            return d
    raise ValueError(f"r = {r} has no order dividing {M} mod {e}")


def phi(n):
    result, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            result -= result // p
        p += 1
    if m > 1:
        result -= result // m
    return result


def layer_classes(N, M, r):
    """per layer b: (d_b, #classes in layer b) by the divisor formula sum_{e | d_b} phi(e) / ord_e(r)."""
    divs = divisors(N)
    phi_over_ord = {e: (phi(e), order_mod(r, e, M)) for e in divs}
    out = []
    for b in range(M):
        d = math.gcd(N, (pow(r, b, N) - 1) % N) if N > 1 else 1
        if d == 0:
            d = N
        cnt = 0
        for e in divs:
            if d % e == 0:
                ph, o = phi_over_ord[e]
                assert ph % o == 0, (N, M, r, e, ph, o)
                cnt += ph // o
        out.append((d, cnt))
    return out


def class_count(N, M, r):
    """k(G) for G = Z_N x|_r Z_M (requires r^M = 1 mod N): the layer sum."""
    assert pow(r, M, N) == 1 % N, f"r^M != 1 mod N for N={N} M={M} r={r}"
    return sum(c for _, c in layer_classes(N, M, r))


def class_count_closed(N, M, r, divs=None):
    """k(G) = M * sum_{e | N} phi(e) / ord_e(r)^2 (exact in integers: M phi(e) is divisible by ord_e(r)^2 termwise? not
    necessarily, so the sum is formed over a common denominator).  `divs` may pass the divisors of N when known."""
    from fractions import Fraction
    assert pow(r, M, N) == 1 % N, f"r^M != 1 mod N for N={N} M={M} r={r}"
    total = Fraction(0)
    for e in (divs or divisors(N)):
        o = order_mod(r, e, M)
        total += Fraction(phi(e), o * o)
    k = total * M
    assert k.denominator == 1, (N, M, r, k)
    return int(k)


def class_terms(N, M, r):
    """the terms of the closed form: [(e, phi(e), ord_e(r), M phi(e) / ord_e(r)^2)] for e | N, largest contribution first."""
    out = []
    for e in divisors(N):
        o = order_mod(r, e, M)
        out.append((e, phi(e), o, M * phi(e) / (o * o)))
    return sorted(out, key=lambda t: -t[3])


def class_count_burnside(N, M, r):
    """the same by Burnside per layer: (1/o_b) sum_k gcd(d_b, r^k - 1)."""
    total = 0
    for b in range(M):
        d = math.gcd(N, (pow(r, b, N) - 1) % N) if N > 1 else 1
        if d == 0:
            d = N
        o = order_mod(r, d, M)
        s = sum(math.gcd(d, (pow(r, k, d) - 1) % d) if d > 1 else 1 for k in range(o))
        assert s % o == 0
        total += s // o
    return total


def commuting_probability(N, M, r):
    return class_count(N, M, r) / (N * M)


def brute_force_classes(N, M, r):
    """orbit enumeration of the conjugation action on all N M elements (a, b) <-> x^a y^b."""
    n = N * M
    rp = [pow(r, e, N) for e in range(M)]

    def mul(g, h):
        a1, b1 = g
        a2, b2 = h
        return ((a1 + rp[b1] * a2) % N, (b1 + b2) % M)

    def inv(g):
        a, b = g
        bi = (-b) % M
        return ((-rp[bi] * a) % N, bi)

    elems = [(a, b) for b in range(M) for a in range(N)]
    seen = set()
    classes = 0
    for g in elems:
        if g in seen:
            continue
        classes += 1
        orbit = set()
        for h in elems:
            orbit.add(mul(mul(h, g), inv(h)))
        seen |= orbit
    assert len(seen) == n
    return classes


def brute_force_commuting_pairs(N, M, r):
    """#{(g, h) : g h = h g} / |G|^2, an independent route to cp (k(G) = #commuting pairs / |G|)."""
    rp = [pow(r, e, N) for e in range(M)]
    elems = [(a, b) for b in range(M) for a in range(N)]
    c = 0
    for (a1, b1) in elems:
        for (a2, b2) in elems:
            if (a1 + rp[b1] * a2) % N == (a2 + rp[b2] * a1) % N:
                c += 1
    return c


def selftest():
    cases = []
    # (N, M, r): even M (faithful: ord r = M), composite N, non-faithful, odd M, dihedral, abelian
    cases += [(7, 6, 3), (7, 6, 5), (9, 6, 2), (15, 4, 2), (15, 4, 7), (21, 6, 2), (21, 6, 5), (21, 6, 10), (35, 12, 3), (39, 12, 2),
              (45, 12, 2), (55, 20, 2), (33, 10, 2), (65, 12, 2), (91, 12, 3), (105, 12, 2), (13, 12, 2), (13, 6, 4), (13, 4, 5),
              (15, 8, 2), (51, 16, 5), (85, 16, 3), (5, 2, 4), (9, 2, 8), (15, 2, 14), (27, 6, 2), (25, 20, 2), (25, 10, 4),
              (11, 5, 3), (31, 5, 2), (7, 3, 2), (21, 3, 4), (7, 1, 1), (12, 4, 1), (16, 4, 3), (16, 2, 7), (20, 4, 3), (8, 2, 3),
              (8, 2, 5), (9, 3, 4), (63, 6, 2), (63, 6, 5), (65, 4, 8), (35, 12, 2), (33, 10, 10), (15, 4, 4), (21, 6, 8)]
    bad = 0
    tested = 0
    for N, M, r in cases:
        if pow(r, M, N) != 1 % N:
            print(f"skip N={N} M={M} r={r}: r^M != 1")
            continue
        tested += 1
        k1 = class_count(N, M, r)
        k2 = class_count_burnside(N, M, r)
        k3 = class_count_closed(N, M, r)
        kb = brute_force_classes(N, M, r)
        cp_pairs = brute_force_commuting_pairs(N, M, r)
        ordr = order_mod(r, N, M)
        faithful = "faithful" if ordr == M else f"ord r = {ordr}"
        inv = sum(1 for b in (0, M // 2) if M % 2 == 0 or b == 0 for a in range(N)
                  if (a * (1 + pow(r, b, N))) % N == 0 and (2 * b) % M == 0 and (a, b) != (0, 0)) if M % 2 == 0 else 0
        ok = (k1 == k2 == k3 == kb) and (cp_pairs == kb * N * M)
        bad += not ok
        print(f"N={N:3d} M={M:2d} r={r:2d} |G|={N * M:5d} {faithful:12s} involutions={inv:3d}  k(G): layer-sum {k1:4d}  burnside {k2:4d}  "
              f"closed-form {k3:4d}  brute-force {kb:4d}  commuting pairs/|G| {cp_pairs // (N * M):4d}  cp={kb / (N * M):.4f}  {'OK' if ok else 'MISMATCH'}")
    # closed-form checks
    for N in (7, 11, 13, 23):
        for M in divisors(N - 1):
            if M < 2:
                continue
            # a primitive root power of order M
            g = next(g for g in range(2, N) if order_mod(g, N, N - 1) == N - 1)
            r = pow(g, (N - 1) // M, N)
            k = class_count(N, M, r)
            assert k == (N - 1) // M + M, (N, M, r, k)
    for N in (5, 7, 9, 15, 21):
        assert class_count(N, 2, N - 1) == (N + 3) // 2
        assert class_count(N, 4, 1) == 4 * N
    # class counts computed independently on the full multiplication table of four groups
    table_counts = {(1275, 16, 107): 171, (1025, 20, 98): 106, (341, 60, 18): 94, (623, 33, 2): 95}
    for (N, M, r), k in table_counts.items():
        got = class_count(N, M, r)
        got2 = class_count_closed(N, M, r)
        print(f"table cross-check Z_{N} x|_{r} Z_{M}: k(G) = {got} / closed form {got2} (full table: {k}) {'OK' if got == got2 == k else 'MISMATCH'}")
        bad += not (got == got2 == k)
    # the layer-pair form cp = sum_{b,d} gcd(d_b, d_d) / (N M^2) against k(G) / |G|
    for N, M, r in ((135269, 72, 8594), (21793, 90, 34), (1275, 16, 107), (45, 12, 2), (105, 12, 2), (15, 4, 7), (33, 10, 10)):
        dd = [math.gcd(N, (pow(r, b, N) - 1) % N) or N for b in range(M)]
        s_pairs = sum(math.gcd(dd[b], dd[e]) for b in range(M) for e in range(M))
        ok = s_pairs == class_count_closed(N, M, r) * M   # cp N M^2 = k(G) M
        print(f"layer-pair form Z_{N} x|_{r} Z_{M}: sum gcd(d_b, d_d) = {s_pairs} = M k(G) = {class_count_closed(N, M, r) * M} {'OK' if ok else 'MISMATCH'}")
        bad += not ok
    # closed form against the layer sum on larger groups (no brute force there)
    for N, M, r in ((135269, 72, 8594), (512257, 160, 706), (633313, 128, 3996), (24461, 150, 16), (86731, 112, 3094), (21793, 90, 34)):
        a, b = class_count(N, M, r), class_count_closed(N, M, r)
        print(f"large-group check Z_{N} x|_{r} Z_{M}: layer sum {a}, closed form {b} {'OK' if a == b else 'MISMATCH'}")
        bad += a != b
    print(f"{tested - bad} of {tested} brute-force cases agree (formula = Burnside = brute-force orbits = commuting pairs / |G|); "
          f"closed forms (prime N, dihedral, abelian) OK" if bad == 0 else f"{bad} MISMATCHES")
    return bad == 0


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        sys.exit(0 if selftest() else 1)
    N, M, r = (int(v) for v in sys.argv[1:4])
    layers = layer_classes(N, M, r)
    k = sum(c for _, c in layers)
    print(f"G = Z_{N} x|_{r} Z_{M}, |G| = {N * M}, ord r = {order_mod(r, N, M)}, k(G) = {k}, cp = {k / (N * M):.6g}")
    print("layer b: d_b = gcd(N, r^b - 1), classes:")
    for b, (d, c) in enumerate(layers):
        if d > 1 or b == 0:
            print(f"  b={b:4d} d_b={d:8d} classes={c}")
    print(f"  ({sum(1 for d, _ in layers if d == 1)} layers with d_b = 1: one class each)")
