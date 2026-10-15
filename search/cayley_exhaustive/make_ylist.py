"""Complete level-pattern list for lift's `Y la lb lc` sweeps (3 inverse pairs, Cayley graph on Z_N x|_r Z_M, N prime).

usage: make_ylist.py N M RANK_FILE OUT_PREFIX CHUNK_CANDIDATES      (RANK_FILE = "-" for no ranking)

Enumerates every level multiset la <= lb <= lc with levels in 0..(M-1)/2 (one representative per inverse pair),
at most one level 0 and gcd(la, lb, lc, M) = 1.  Patterns with two generators in Z_N are excluded because the
quotient Cay(Z_M, {+-lc}) is a cycle of diameter floor(M/2) > D for the diameters this is used for.
Order: patterns missing from the ranking first (so that they are visibly covered), then the ranked patterns in the
ranking's order ("PS coinc P : la lb lc" lines as printed by cayley_sa --patscreen: radius-2 coincidences of the Z_M
quotient ascending, then the Poisson yardstick P ascending), then the level-0 patterns; without a ranking file every
pattern with la > 0 comes first in lexicographic order.
A pattern with la > 0 costs N + 2 candidates in lift (a = y^la; b = x y^lb, c = x^t y^lc for all t; b = y^lb,
c in {y^lc, x y^lc}); a pattern with la = 0 costs N (a = x, b = x y^lb, c = x^t y^lc for all t).
Writes OUT_PREFIX_c<k>.txt chunk files of about CHUNK_CANDIDATES candidates each and OUT_PREFIX.manifest with the
exact pattern and candidate counts per chunk (a completeness certificate compares lift's EXHAUSTIVE totals against it).
"""
import math
import sys


def main():
    N, M, rankfn, outp, chunk = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4], int(sys.argv[5])
    L = (M - 1) // 2
    allpats = set()
    for la in range(0, L + 1):
        for lb in range(max(la, 1), L + 1):
            for lc in range(lb, L + 1):
                if math.gcd(math.gcd(la, lb), math.gcd(lc, M)) == 1:
                    allpats.add((la, lb, lc))
    ranked = []
    seen = set()
    for line in (open(rankfn) if rankfn != "-" else []):
        if not line.startswith("PS "):
            continue
        parts = line.split(":")
        lv = tuple(sorted(int(x) for x in parts[1].split()))
        head = parts[0].split()
        if lv in seen:
            raise SystemExit(f"duplicate pattern in ranking: {lv}")
        if lv not in allpats:
            raise SystemExit(f"ranking has a pattern outside the enumeration: {lv}")
        seen.add(lv)
        ranked.append((lv, int(float(head[1])), float(head[2])))
    nonzero = sorted(p for p in allpats if p[0] > 0)
    zero = sorted(p for p in allpats if p[0] == 0)
    missing = [p for p in nonzero if p not in seen]
    missing_txt = str(missing) if rankfn != "-" else "(no ranking given)"
    order = [(p, None, None) for p in missing] + ranked + [(p, None, None) for p in zero]
    assert len(order) == len(allpats) and len(set(p for p, _, _ in order)) == len(allpats)
    cost = lambda p: N + 2 if p[0] > 0 else N
    total = sum(cost(p) for p in allpats)
    chunks, cur, curc = [], [], 0
    for item in order:
        cur.append(item)
        curc += cost(item[0])
        if curc >= chunk:
            chunks.append((cur, curc))
            cur, curc = [], 0
    if cur:
        chunks.append((cur, curc))
    with open(outp + ".manifest", "w") as mf:
        mf.write(f"# N={N} M={M} L={L}: {len(allpats)} patterns = {len(nonzero)} with la>0 (cost N+2={N + 2} each) "
                 f"+ {len(zero)} with la=0 (cost N={N} each); total candidates {total}; "
                 f"{len(ranked)} ranked, {len(missing)} missing from the ranking: {missing_txt}\n")
        for k, (items, c) in enumerate(chunks):
            fn = f"{outp}_c{k}.txt"
            nz = sum(1 for p, _, _ in items if p[0] > 0)
            with open(fn, "w") as f:
                f.write(f"# chunk {k}: {len(items)} patterns ({nz} with la>0, {len(items) - nz} with la=0), {c} candidates\n")
                for p, coinc, P in items:
                    tail = f"   # coinc {coinc} P {P:.3f}" if coinc is not None else ("   # unranked" if p[0] > 0 else "")
                    f.write(f"Y {p[0]} {p[1]} {p[2]}{tail}\n")
            mf.write(f"{fn} patterns={len(items)} la0={len(items) - nz} candidates={c}\n")
    print(f"{len(allpats)} patterns ({len(nonzero)} + {len(zero)} level-0), {total} candidates, "
          f"{len(chunks)} chunks; missing from ranking: {missing_txt}")


if __name__ == "__main__":
    main()
