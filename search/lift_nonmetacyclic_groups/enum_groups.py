#!/usr/bin/env python3
"""enum_groups.py: the groups G = K x| C_M of the block family (enum_family.py's blocks and automorphism classes) in a
window (LO, HI] with filters for two uses:

  --cell cay7   (7,6) one-involution Cayley graphs: M even, alpha faithful or with kernel exactly <y^(M/2)>, |Z(G) ∩ K| = 1
                (|Z(G)| <= 2), any number of involutions (all involutions are K-conjugate when |K| is odd, so the involution
                of the connection set is y^(M/2) without loss of generality).
  --cell lift   shape-103 two-vertex lifts: any M >= 1 (M = 1: G = K itself, alpha = id), alpha faithful or kernel <y^(M/2)>
                (M even), |Z(K)^alpha| <= ZMAX (default 9); no involution filter.

Output: specs JSONL (id, n, M, K, variant, metacyclic flag, block decomposition, spec) after enum_family's power-map
reduction alpha ~ alpha^j.  Metacyclic specs (all blocks cyclic) are kept and flagged (they are the controls / the
w8-w14 / w20 coverage); the non-metacyclic ones are the new family.

usage: enum_groups.py LO HI --cell cay7|lift [--mmin 1] [--mmax 150] [--zmax 9] --out specs.jsonl
"""
import itertools
import json
import sys
import time
from collections import defaultdict
from math import gcd

import enum_family as ef  # noqa: E402


def main():
    LO, HI = int(sys.argv[1]), int(sys.argv[2])
    arg = lambda k, d: type(d)(sys.argv[sys.argv.index(k) + 1]) if k in sys.argv else d
    cell = arg("--cell", "lift")
    mmin, mmax, zmax = arg("--mmin", 1), arg("--mmax", 150), arg("--zmax", 9)
    outfn = arg("--out", "")
    prefix = arg("--prefix", "grp")
    if cell == "cay7":
        mmin = max(mmin, 2)
    cand = []
    stats = defaultdict(int)
    t0 = time.time()
    for n in range(LO + 1, HI + 1):
        for M in range(mmin, mmax + 1):
            if n % M:
                continue
            if cell == "cay7" and M % 2:
                continue
            K = n // M
            if K % 2 == 0 or K < 3:
                continue
            for dec in ef.decompositions(K):
                for variant in ("faithful", "kernel2"):
                    target = M if variant == "faithful" else (M // 2 if M % 2 == 0 else 0)
                    if target == 0:
                        continue
                    fl = [[a for a in (ef.autos_cyc(bl[1], target) if bl[0] == "cyc" else ef.get_autos(bl)) if target % a["order"] == 0] for bl in dec]
                    if any(len(L) == 0 for L in fl):
                        continue
                    for combo in itertools.product(*fl):
                        o = 1
                        for a in combo:
                            o = ef.lcm(o, a["order"])
                        if o != target:
                            continue
                        stats["combos"] += 1
                        if M % 2 == 0:
                            if variant == "faithful":
                                fixp = 1
                                for a in combo:
                                    fixp *= a["fix"][(M // 2) % a["order"]]
                                ninv = K // fixp
                            else:
                                ninv = 1
                        else:
                            ninv = 0
                        # centre of G inside K = Z(K)^alpha (alpha-fixed centre of K); Z(G) = that x (kernel of the action on y)
                        zf = 1
                        for a in combo:
                            zf *= a["zfix"][1 % a["order"]] if a["order"] > 1 else a["zfix"][0]
                        zG = zf * (M // target)
                        if cell == "cay7" and zf > 1:
                            stats["z_cut"] += 1; continue
                        if cell == "lift" and zG > zmax:
                            stats["z_cut"] += 1; continue
                        cand.append({"n": n, "M": M, "K": K, "variant": variant, "ninv": ninv, "zG": zG, "dec": dec,
                                     "blocks": [a["spec"] for a in combo], "orders": [a["order"] for a in combo],
                                     "twist": any(a["extra"].get("twist") for a in combo),
                                     "metacyclic": all(bl[0] == "cyc" for bl in dec)})
        if n % 25 == 0:
            print(f"# n = {n}: {len(cand)} candidates, {time.time() - t0:.0f}s", file=sys.stderr, flush=True)
    print(f"# combos {stats['combos']}, centre cut {stats['z_cut']}, surviving {len(cand)}", file=sys.stderr)
    kept = []
    seen = set()
    for c in cand:
        if c["twist"]:
            kept.append(c); continue
        M = c["M"]
        forms = []
        for j in range(1, M + 1):
            if gcd(j, M) != 1:
                continue
            forms.append(tuple(ef.pow_key(bl, sp, o, j) for bl, sp, o in zip(c["dec"], c["blocks"], c["orders"])))
        cf = (c["n"], c["M"], c["variant"], tuple(c["dec"]), min(forms))
        if cf in seen:
            continue
        seen.add(cf)
        kept.append(c)
    nm = sum(1 for c in kept if not c["metacyclic"])
    print(f"# after power-map reduction: {len(kept)} ({nm} non-metacyclic, {len(kept) - nm} metacyclic)", file=sys.stderr)
    if outfn:
        with open(outfn, "w") as f:
            for i, c in enumerate(kept):
                spec = {"M": c["M"], "blocks": c["blocks"]}
                f.write(json.dumps({"id": f"{prefix}_{c['n']}_{i}", "n": c["n"], "M": c["M"], "K": c["K"], "variant": c["variant"],
                                    "ninv_pred": c["ninv"], "zG_pred": c["zG"], "metacyclic": c["metacyclic"],
                                    "dec": [list(d) for d in c["dec"]], "twist": c["twist"], "spec": spec}) + "\n")


if __name__ == "__main__":
    main()
