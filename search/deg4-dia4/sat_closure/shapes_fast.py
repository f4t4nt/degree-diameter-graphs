"""Shape enumeration for three K4,4 blocks with slack, n = 72 + N, fast version.

Same mathematics as shapes107.py / shapes_stream.py; the search is reorganised:

  1. components107.enumerate_components streams every closed typed component; a gadget (non-unit
     component) is stored compactly under its fine signature
         (deficit, pot4, e, n, ex)
     with e = slots - n (1,1,1) (sum_j e_j = -deficit) and ex = the 6-vector of per-ordered-pair
     excesses rows_ij - 2 slots_i (sum = 4 potential + deficit... in quarter units: see below).
     Balance: the gadget connectors hold g = G - 3 ports per block (G gadget connectors, deficit
     D = 3N - 96 = 9 at 107), i.e. sum over gadgets of e = -3 (1,1,1); rows: sum of ex >= 0 in
     every coordinate (units contribute exactly 2 per slot).
  2. Coarse recursion over keys (deficit, pot4, e, n) (multisets with replacement, deficit sum D
     plus quads, potential >= bound, sum e = (-3,-3,-3), N - conn even), pruned per coordinate by
     a table of the minimum and maximum e_j still reachable with the remaining deficit and cost.
  3. For every coarse multiset, a fine sub-search chooses an ex-vector per element with total >= 0
     (pruned by the Pareto maxima of the remaining choices); then expansion to gadget multisets,
     the port-lemma matching (inequality; equality for tight shapes) from shapes107.py.
  4. Shapes are deduplicated under the S3 of block relabellings (certificates of the relabelled
     gadgets), which the SAT model does not need to distinguish.

Usage: shapes_fast.py N --out=shapes.json [--log=FILE]
"""
import collections
import itertools
import json
import sys
import time
from fractions import Fraction as FR

import components107 as C
import shapes107 as S

PAIRS = S.PAIRS
PERMS = list(itertools.permutations(range(3)))


def analyse_raw(types, edges):
    n = len(types)
    nbr = {i: set() for i in range(n)}
    for a, b in edges:
        nbr[a].add(b)
        nbr[b].add(a)
    f = [tuple(types[i][j] + sum(types[k][j] for k in nbr[i]) for j in range(3)) for i in range(n)]
    deficit = sum(3 - sum(t) for t in types)
    pot = sum(C.ETAB[(types[a], types[b])] for a, b in edges) if edges else C.iso_pot(types[0])
    pot4 = int(pot * 4)
    slots = tuple(sum(t[j] for t in types) for j in range(3))
    rows = tuple(sum(types[h][i] * f[h][j] for h in range(n)) for (i, j) in PAIRS)
    ex = tuple(rows[k] - 2 * slots[PAIRS[k][0]] for k in range(6))
    e = tuple(slots[j] - n for j in range(3))
    return (deficit, pot4, e, n, ex)


def pack(types, edges):
    return bytes(v for t in types for v in t), bytes(v for e in sorted(edges) for v in e)


def unpack(tb, eb):
    types = [tuple(tb[3 * i:3 * i + 3]) for i in range(len(tb) // 3)]
    edges = [(eb[2 * i], eb[2 * i + 1]) for i in range(len(eb) // 2)]
    return types, edges


def gadget_record(tb, eb):
    types, edges = unpack(tb, eb)
    c = {"types": [list(t) for t in types], "edges": [list(e) for e in edges], "n": len(types),
         "slots": [sum(t[j] for t in types) for j in range(3)]}
    pot = sum(C.ETAB[(types[a], types[b])] for a, b in edges) if edges else C.iso_pot(types[0])
    c["potential"] = str(pot)
    return S.analyse(c)


def shape_key(gadgets):
    """Canonical key under block relabelling: min over S3 of the sorted gadget certificates."""
    best = None
    for p in PERMS:
        certs = []
        for g in gadgets:
            types = [tuple(t[p[j]] for j in range(3)) for t in g["types"]]
            certs.append(C.closed_cert(types, {tuple(e) for e in g["edges"]}))
        k = tuple(sorted(certs))
        if best is None or k < best:
            best = k
    return best


def main():
    t0 = time.time()
    N = int(sys.argv[1])
    out = next((a.split("=")[1] for a in sys.argv if a.startswith("--out=")), None)
    logp = next((a.split("=")[1] for a in sys.argv if a.startswith("--log=")), None)
    log = open(logp, "a") if logp else sys.stderr
    family = next((a.split("=")[1] for a in sys.argv if a.startswith("--family=")), None)
    family = int(family) if family is not None else None
    bound4 = (96 - 3 * N) * 4 if family is None else 4 * family
    Dtot = 3 * N - 96
    store = collections.defaultdict(list)       # fine signature -> packed gadgets
    stats = collections.Counter()

    def sink(types, edges, pot):
        if S.is_unit({"n": len(types), "types": types}):
            stats["units"] += 1
            return
        sig = analyse_raw(types, edges)
        if sig[1] < bound4:
            stats["below_bound"] += 1
            return
        store[sig].append(pack(types, edges))
        stats["gadgets"] += 1

    summary = C.enumerate_components(N, log=log, sink=sink, family_bound=family)
    print(json.dumps(summary), file=log, flush=True)
    coarse = collections.defaultdict(list)      # (deficit, pot4, e, n) -> list of fine signatures
    for sig in store:
        coarse[sig[:4]].append(sig)
    keys = sorted(coarse, key=lambda k: (k[0], -k[1]))
    print(f"{stats['gadgets']} gadgets, {len(store)} fine signatures, {len(keys)} coarse keys, {time.time() - t0:.0f} s",
          file=log, flush=True)
    maxdef = summary["maxdef"]
    # Range table: for multisets of coarse keys with deficit sum exactly r (quads count -1, so
    # r >= -2) and cost exactly p (quarter units), the componentwise minimum and maximum of the
    # summed e vector; unbounded knapsack relaxation (components107's bound keeps every item cost
    # >= 0 and the table finite: deficit-0 items cost >= 8, zero-cost items have deficit >= 4).
    Pmax = -bound4
    INF = 10 ** 9
    items = [(k[0], -k[1], k[2]) for k in keys]
    lo3 = {(r, p): [INF] * 3 for r in range(-2, maxdef + 1) for p in range(Pmax + 1)}
    hi3 = {(r, p): [-INF] * 3 for r in range(-2, maxdef + 1) for p in range(Pmax + 1)}
    lo3[(0, 0)] = [0, 0, 0]
    hi3[(0, 0)] = [0, 0, 0]
    changed = True
    while changed:
        changed = False
        for (r, p) in sorted(lo3, key=lambda k: (k[1], k[0])):
            if lo3[(r, p)][0] == INF:
                continue
            for d, c, e in items:
                r2, p2 = r + d, p + c
                if (r2, p2) not in lo3:
                    continue
                for j in range(3):
                    v = lo3[(r, p)][j] + e[j]
                    if v < lo3[(r2, p2)][j]:
                        lo3[(r2, p2)][j] = v
                        changed = True
                    v = hi3[(r, p)][j] + e[j]
                    if v > hi3[(r2, p2)][j]:
                        hi3[(r2, p2)][j] = v
                        changed = True
    print(f"range table done, {time.time() - t0:.0f} s", file=log, flush=True)

    def feasible_rest(R, P, need):
        """Can a multiset with deficit sum exactly R and cost <= P have e-sum equal to need?
        (Necessary condition: need within the componentwise range for some cost p <= P.)"""
        for p in range(P + 1):
            k = (R, p)
            if k in lo3 and lo3[k][0] != INF and all(lo3[k][j] <= need[j] <= hi3[k][j] for j in range(3)):
                return True
        return False

    last_quad = max([i for i, k in enumerate(keys) if k[0] < 0], default=-1)
    coarse_shapes = []
    target = (-Dtot // 3,) * 3  # sum of e over gadgets: g - G = -D/3
    nosym = "--nosym" in sys.argv

    def rec(start, chosen, pot4, deff, e, conn):
        stats["rec_nodes"] += 1
        if deff == Dtot and chosen and e == target:
            rest = N - conn
            if rest % 2 == 0:
                coarse_shapes.append((rest // 2, pot4, list(chosen)))
                stats["coarse_balanced"] += 1
        P = pot4 - bound4
        quads_ok = (P // 20) if start <= last_quad else 0
        for k in range(start, len(keys)):
            key = keys[k]
            d, p4, ek, n = key
            if pot4 + p4 < bound4:
                continue
            nd = deff + d
            if nd > Dtot + quads_ok:
                continue
            if conn + n > N:
                continue
            ne = (e[0] + ek[0], e[1] + ek[1], e[2] + ek[2])
            need = (target[0] - ne[0], target[1] - ne[1], target[2] - ne[2])
            if not feasible_rest(Dtot - nd, pot4 + p4 - bound4, need):
                stats["pruned_range"] += 1
                continue
            rec(k, chosen + [key], pot4 + p4, nd, ne, conn + n)

    rec(0, [], 0, 0, (0, 0, 0), 0)
    print(json.dumps(dict(stats)), f"{time.time() - t0:.0f} s", file=log, flush=True)

    # fine sub-search: choose ex-vectors with total >= 0
    fine_shapes = []
    for u, pot4, chosen in coarse_shapes:
        cnt = collections.Counter(chosen)
        slots = []  # list of (key, multiplicity)
        for key, m in cnt.items():
            slots.append((key, m))
        options = [[(sig[4], sig) for sig in coarse[key]] for key, m in slots]
        # Pareto maxima per slot for pruning: componentwise max over options times multiplicity
        upper = [tuple(max(o[0][j] for o in opts) * m for j in range(6)) for opts, (key, m) in zip(options, slots)]
        suffix = [(0,) * 6] * (len(slots) + 1)
        for i in range(len(slots) - 1, -1, -1):
            suffix[i] = tuple(suffix[i + 1][j] + upper[i][j] for j in range(6))

        def frec(i, total, picked):
            if i == len(slots):
                if all(v >= 0 for v in total):
                    fine_shapes.append((u, pot4, list(picked)))
                    stats["fine_ok"] += 1
                return
            if any(total[j] + suffix[i][j] < 0 for j in range(6)):
                stats["fine_pruned"] += 1
                return
            key, m = slots[i]
            for combo in itertools.combinations_with_replacement(options[i], m):
                t2 = list(total)
                for ex, sig in combo:
                    for j in range(6):
                        t2[j] += ex[j]
                frec(i + 1, tuple(t2), picked + [sig for _, sig in combo])

        frec(0, (0,) * 6, [])
    print(json.dumps(dict(stats)), f"{time.time() - t0:.0f} s", file=log, flush=True)

    shapes = []
    cache = {}
    seen_keys = set()
    for u, pot4, sigs in fine_shapes:
        cnt = collections.Counter(sigs)
        for s in cnt:
            if s not in cache:
                cache[s] = [gadget_record(tb, eb) for tb, eb in store[s]]
        pools = [list(itertools.combinations_with_replacement(cache[s], m)) for s, m in cnt.items()]
        for combo in itertools.product(*pools):
            gs = [g for part in combo for g in part]
            stats["gadget_multisets"] += 1
            if not S.matching_ok(gs):
                stats["matching_fail"] += 1
                continue
            if pot4 == 4 * (96 - 3 * N) and not S.matching_ok(gs, exact=True):
                stats["tight_exact_matching_fail"] += 1
                continue
            if not nosym:
                k = shape_key(gs)
                if k in seen_keys:
                    stats["block_relabel_dup"] += 1
                    continue
                seen_keys.add(k)
            stats["survive"] += 1
            shapes.append({"units": u, "potential": str(FR(pot4, 4)), "budget": 768 + (pot4 - 4 * (96 - 3 * N)) // 2,
                           "gadgets": [{"types": g["types"], "edges": g["edges"], "potential": g["potential"], "f": g["f"]} for g in gs]})
    print(json.dumps(dict(stats)), f"{time.time() - t0:.0f} s", file=log, flush=True)
    shapes.sort(key=lambda s: (-FR(s["potential"]), len(s["gadgets"]), sum(len(g["types"]) for g in s["gadgets"])))
    result = {"N": N, "n": 72 + N, "bound": str(FR(bound4, 4)), "enumeration": summary, "stats": dict(stats),
              "fine_signatures": len(store), "coarse_keys": len(keys), "wall_s": round(time.time() - t0, 1), "shapes": shapes}
    if out:
        json.dump(result, open(out, "w"), indent=0)
    else:
        json.dump(result, sys.stdout, indent=0)


if __name__ == "__main__":
    main()
