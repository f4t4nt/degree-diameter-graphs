"""Brute-force completeness check of gen_cay.py on a small group Z_N x|_r Z_M (odd order, or even M: involutions are
skipped as generators since they never give a simple 6-regular Cayley graph with three inverse pairs, and lift.c
rejects them in is_simple()).

Every inverse-closed 6-set S = {a^±1, b^±1, c^±1} is canonicalised under the explicit automorphism group
A = {(x,y) -> (i x + j S_y, y)} (all units i, all valid j); the set of canonical forms met by the generated
candidates (E lines expanded exactly as lift.c does) must equal the set of all canonical forms.

usage: check_cay.py N M r
"""
import subprocess
import sys
from itertools import combinations
from math import gcd

HERE = __file__.rsplit("/", 1)[0] if "/" in __file__ else "."


def main():
    N, M, r = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    m = N * M
    H = (M - 1) // 2
    rp = [pow(r, y, N) for y in range(M)]
    Sy = [0] * (M + 1)
    for y in range(M):
        Sy[y + 1] = (Sy[y] + rp[y]) % N
    j0 = N // gcd(N, Sy[M])
    J = list(range(0, N, j0))
    units = [i for i in range(1, N) if gcd(i, N) == 1]

    def inv(g):
        x, y = g % N, g // N
        yi = (M - y) % M
        return (-rp[yi] * x) % N + N * yi

    def mul(g, h):
        x1, y1, x2, y2 = g % N, g // N, h % N, h // N
        return (x1 + rp[y1] * x2) % N + N * ((y1 + y2) % M)

    # explicit automorphisms as permutations; verify the homomorphism property on a sample
    auts = []
    for i in units:
        for j in J:
            perm = [((i * (g % N) + j * Sy[g // N]) % N) + N * (g // N) for g in range(m)]
            auts.append(perm)
    import random

    random.seed(1)
    for perm in random.sample(auts, min(20, len(auts))):
        for _ in range(50):
            g, h = random.randrange(m), random.randrange(m)
            assert perm[mul(g, h)] == mul(perm[g], perm[h]), "not an automorphism"
        assert len(set(perm)) == m

    def canon_set(S):
        return min(tuple(sorted(perm[s] for s in S)) for perm in auts)

    pairs = []
    seen = set()
    for g in range(1, m):
        if g in seen:
            continue
        gi = inv(g)
        seen |= {g, gi}
        if gi == g:
            continue  # involution: cannot be one of three inverse pairs of a simple 6-regular Cayley graph
        pairs.append((g, gi))
    all_forms = set()
    n_disc = 0
    for P in combinations(pairs, 3):
        S = [s for pr in P for s in pr]
        # sets whose y-coordinates do not generate Z_M give disconnected Cayley graphs and are skipped by design
        g = M
        for s in S:
            g = gcd(g, s // N)
        if g != 1:
            n_disc += 1
            continue
        all_forms.add(canon_set(S))
    print(f"6-sets with y-coordinates in a proper subgroup of Z_M (skipped by design): {n_disc}")

    res = subprocess.run([sys.executable, f"{HERE}/gen_cay.py", str(N), str(M), str(r)], capture_output=True, text=True, check=True)
    cands = []
    for line in res.stdout.splitlines():
        f = line.split()
        if f[0] == "T":
            cands.append((int(f[1]), int(f[2]), int(f[3])))
        elif f[0] == "E":
            a, b, L = int(f[1]), int(f[2]), int(f[3])
            gy = gcd(gcd(a // N, b // N), M)
            if L == 0:
                if gy == 1:
                    cands.extend((a, b, x) for x in range(1, (N - 1) // 2 + 1))
                L = 1
            for y in range(L, H + 1):
                if gcd(gy, y) != 1:
                    continue
                cands.extend((a, b, x + N * y) for x in range(N))
            # mirrors lift.c enumerate_last: level M/2 for even M, canonical under inversion
            if M % 2 == 0 and L <= M // 2 and gcd(gy, M // 2) == 1:
                y = M // 2
                cands.extend((a, b, x + N * y) for x in range(N) if inv(x + N * y) > x + N * y)
    hit = set()
    for a, b, c in cands:
        S = {a, inv(a), b, inv(b), c, inv(c)}
        if len(S) != 6:
            continue
        hit.add(canon_set(S))
    missing = all_forms - hit
    print(f"N={N} M={M} r={r} m={m} |A|={len(auts)} 6-sets={sum(1 for _ in combinations(pairs, 3))} orbits={len(all_forms)} "
          f"candidates={len(cands)} orbits hit={len(hit)} missing={len(missing)} over-coverage={len(cands) / max(1, len(all_forms)):.2f}")
    print(res.stderr.strip())
    if missing:
        print("MISSING, e.g.", sorted(missing)[:5])
        sys.exit(1)


if __name__ == "__main__":
    main()
