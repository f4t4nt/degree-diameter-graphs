"""Demand-driven window-multiset enumeration (clusters) for three K4,4 blocks, n = 72 + N.

Windows, costs and shape quantities as in windows108.py (imported). A free DFS over windows is
hopeless (2M nodes without finishing at 106), so the enumeration grows multisets from a seed
window by *demands*:

  - half-edge demand: a(t -> m) > a(m -> t) (or a(t -> t) odd) means the final multiset contains
    further windows of type m carrying t (the chaining condition is an equality in every member);
    these demands are always branched on (all satisfiers of id >= seed);
  - slot demand: eff_j(a) > eff_j(4 - a) for a foreign vector a != (2,2) means a block-j slot of
    vector a still lacks its co-holder. With --kmax=K a slot demand with at most K exact satisfiers
    is branched on (satisfiers, or "declared" inexact when the excess budget allows); demands with
    more satisfiers are left to the combine stage, where the partners may come from other
    clusters. The runs reported in RESULTS.md use --kmax=0: only slots with *no* possible
    complement anywhere are branched (they must be declared inexact, at most 2 X of them per
    cluster, X = excess budget), so a cluster is exactly a chaining-closed window multiset.
    (Values K > 0 are kept for experiments; they are not complete for the decomposition into
    clusters, see RESULTS.md section 2.2.)

Growth from a seed window s (the window of smallest id in the cluster) adds only windows of id >=
s.id and is memoised by state digest; every chaining-closed multiset is reached from its minimum
window. Bounds (all necessary conditions for the whole member): cost + chaining bound (half-edge
partner shares; per demanded type the number of windows needed times their minimum cost), cost +
slot-complement bound (slots without present complement need partners at min cost per slot,
except 2 per unit of unused excess), deficit including affordable quads, demanded-type deficit,
slots <= 32, size <= --maxsize (default N).

Stage combine (plain knapsack with repetition over clusters; fine at 106, too slow above 10^4
clusters; see combine108.py and pairs108.py for the stages used at 108): deficit D (quads -1),
cost in the accepted set, declared slots <= 2 (D - cost), slot balance (32 - 2u) x 3, rows >= 0,
hub-to-hub, (2,2)-parity when exact, imbalance of the slot union <= 2 (D - cost), canonical under
block relabelling, slot-vector matching relaxation; output in the format read by windows108.py
shapes.

Usage:
  clusters108.py clusters N --cmax=C --excess=X [--kmax=K] [--maxsize=S] --out=runs/C.jsonl [--log=FILE] [--closure]
  clusters108.py combine N --clusters=runs/C.jsonl --accept=10,12 --out=runs/W.jsonl [--log=FILE]
"""
import collections
import hashlib
import json
import sys
import time
from fractions import Fraction as FR

import windows108 as W

C = W.C
S = W.S
NONISO = W.NONISO
TIDX = W.TIDX
ETAB = W.ETAB
BIG = W.BIG
NT = len(NONISO)
DEG = [4 - sum(t) for t in NONISO]
DEF = [3 - sum(t) for t in NONISO]
BEST = {t: max(ETAB[(t, y)] for y in NONISO) for t in NONISO}
HC = [[int(-ETAB[(NONISO[t], NONISO[m])] / 2 * 8) for m in range(NT)] for t in range(NT)]


def committed8(w):
    """Sound lower bound on the member's cost attributable to the edges of w and the further edges
    of its neighbours (each further edge may be shared by two neighbours, hence the halves)."""
    if w.quad:
        return w.cost8
    c = 2 * w.cost8
    for m in w.M:
        c += int(-BEST[m] * (4 - sum(m) - 1) * 8) / 2
    return c


def window_table(cmax8, closure, log):
    wins = [w for w in W.all_windows(cmax8) if committed8(w) <= cmax8]
    print(f"windows cost <= {cmax8 / 8}: committed-cost filter leaves {len(wins)}", file=log, flush=True)
    if closure:
        # tight closure: every slot complementable and every half-edge matchable inside the set
        while True:
            prov = set()
            hv = set()
            for w in wins:
                for (j, a) in w.svec:
                    prov.add((j, a))
                if not w.quad:
                    for m in w.M:
                        hv.add((w.x, m))
            keep = []
            for w in wins:
                ok = all(4 - a[0] >= 0 and 4 - a[1] >= 0 and (j, (4 - a[0], 4 - a[1])) in prov for (j, a) in w.svec)
                ok = ok and (w.quad or all((m, w.x) in hv for m in w.M))
                if ok:
                    keep.append(w)
            if len(keep) == len(wins):
                break
            wins = keep
        print(f"tight closure leaves {len(wins)}", file=log, flush=True)
    # re-index
    for k, w in enumerate(wins):
        w.id = k
    return wins


class Growth:
    """Demand-driven growth with explicit slot resolution.

    State = (window counts, resolved slots, excess used). Effective slot counts are the window
    slot counts minus the resolved ones. Demands: half-edge (t -> m) with a(t -> m) > a(m -> t) (or
    a(t -> t) odd), and slot (j, a) with eff_j(a) > eff_j(4 - a), a != (2,2). In a member every
    demand is met in exactly one of these ways, which are the branches:
      half-edge: a further window of type m carrying t;
      slot (j, a): (i) a further window with a block-j slot 4 - a (exact); (ii) a further window
      with a block-j slot b >= 4 - a, b != 4 - a, resolving (a, b) at excess |a + b - 4|; (iii) a
      present unmatched block-j slot b with a + b >= 4 inexactly; (iv) a unit or (2,2) partner when
      a >= (2,2), at excess |a| - 4. The total excess of a member is Phi + D (0 at 768, 2 at 772).
    The demand with the fewest options is branched on. States are memoised by digest."""

    def __init__(self, wins, N, cmax8, X, maxsize, log, kmax=3):
        self.wins = wins
        self.kmax = kmax
        self.N = N
        self.Dtot = 3 * N - 96
        self.cmax8 = cmax8
        self.X = X
        self.maxsize = maxsize
        self.log = log
        self.by_half = collections.defaultdict(list)   # (m, t): windows of type m carrying t
        self.by_slot = collections.defaultdict(list)   # (j, b): windows with a block-j slot of vector b
        for w in wins:
            if not w.quad:
                for m, _ in w.halves:
                    self.by_half[(w.t, m)].append(w)
            for (j, b) in set(w.svec):
                self.by_slot[(j, b)].append(w)
        self.vectors = [sorted({b for (j, b) in self.by_slot if j == jj}) for jj in range(3)]
        self.minsc = {}
        for w in wins:
            if w.nx == 0:
                continue
            c = w.cost8 / w.nx
            for (j, b) in set(w.svec):
                if (j, b) not in self.minsc or c < self.minsc[(j, b)]:
                    self.minsc[(j, b)] = c
        self.minw = [[BIG] * NT for _ in range(NT)]
        self.maxmult = [[0] * NT for _ in range(NT)]
        for w in wins:
            if w.quad:
                continue
            for m, k in w.halves:
                if w.cost8 < self.minw[w.t][m]:
                    self.minw[w.t][m] = w.cost8
                if k > self.maxmult[w.t][m]:
                    self.maxmult[w.t][m] = k
        quad_costs = sorted(w.cost8 for w in wins if w.quad)
        self.cheapest_quad8 = quad_costs[0] if quad_costs else BIG
        self.stats = collections.Counter()
        self.visited = set()

    def partner_cost(self, j, a):
        b = (4 - a[0], 4 - a[1])
        if b[0] < 0 or b[1] < 0:
            return BIG
        return self.minsc.get((j, b), BIG)

    def analyse(self, counts, resolved):
        a = [[0] * NT for _ in range(NT)]
        sv = [collections.Counter() for _ in range(3)]
        slots = [0, 0, 0]
        ex = [0] * 6
        hh = [0, 0, 0]
        st = {"deficit": 0, "cost8": 0, "G": 0, "quads": 0}
        for wid, cnt in counts.items():
            self.add(self.wins[wid], cnt, a, sv, slots, ex, hh, st)
        for (j, b), cnt in resolved.items():
            sv[j][b] -= cnt
        return a, sv, slots, ex, hh, st

    @staticmethod
    def add(w, cnt, a, sv, slots, ex, hh, st):
        if not w.quad:
            for m, k in w.halves:
                a[w.t][m] += k * cnt
        for (j, fv) in w.svec:
            sv[j][fv] += cnt
        for j in range(3):
            slots[j] += cnt * w.x[j]
            hh[j] += cnt * w.hh[j]
        for k in range(6):
            ex[k] += cnt * w.ex[k]
        st["deficit"] += cnt * w.deficit
        st["cost8"] += cnt * w.cost8
        st["G"] += cnt
        if w.quad:
            st["quads"] += cnt

    def bounds_ok(self, a, sv, slots, st, used):
        cost8 = st["cost8"]
        if cost8 > self.cmax8 or st["G"] > self.maxsize or max(slots) > 32:
            return False
        qmax = (self.cmax8 - cost8) // self.cheapest_quad8 if self.cheapest_quad8 < BIG else 0
        if st["deficit"] > self.Dtot + qmax:
            return False
        budget = self.cmax8 - cost8
        # the excess still to be used must be affordable: excess total = D - cost
        if used > self.Dtot - (cost8 + 7) // 8:
            return False
        chain8 = 0
        need_def = 0
        type_bound = 0
        for m in range(NT):
            am = a[m]
            D = am[m] % 2
            nwin = D
            minc = self.minw[m][m] if D else BIG
            if D:
                chain8 += HC[m][m]
            for t in range(NT):
                if t == m:
                    continue
                d = a[t][m] - am[t]
                if d > 0:
                    chain8 += d * HC[t][m]
                    mm = self.maxmult[m][t]
                    if mm == 0:
                        return False
                    nwin = max(nwin, -(-d // mm))
                    if self.minw[m][t] < minc:
                        minc = self.minw[m][t]
                    D += d
            if D:
                nwin = max(nwin, -(-D // DEG[m]))
                type_bound += nwin * minc
                need_def += DEF[m] * nwin
        if chain8 > budget or type_bound > budget:
            return False
        if need_def > self.Dtot - st["deficit"] + qmax:
            return False
        # slot-complement bound: unmatched effective slots need exact partners from future windows
        # except those resolved inexactly later (each unit of excess resolves at most two slots)
        tol = 2 * min(self.X - used, self.Dtot - (cost8 + 7) // 8 - used)
        if tol < 0:
            return False
        costs = []
        for j in range(3):
            c = sv[j]
            for fv, cnt in c.items():
                if fv == (2, 2) or cnt <= 0:
                    continue
                other = c.get((4 - fv[0], 4 - fv[1]), 0)
                if cnt > other:
                    costs.extend([self.partner_cost(j, fv)] * (cnt - other))
        if len(costs) > tol:
            costs.sort(reverse=True)
            if sum(costs[tol:]) > budget + 1e-9:
                return False
        return True

    def grow_all(self, out):
        t0 = time.time()
        fout = open(out, "w")
        for seed in self.wins:
            self.grow(seed, fout)
            self.stats["seeds"] += 1
            if self.stats["seeds"] % 100 == 0:
                print(f"seed {seed.id}/{len(self.wins)}: {json.dumps(dict(self.stats))}, visited {len(self.visited)}, {time.time() - t0:.0f} s",
                      file=self.log, flush=True)
        fout.close()

    @staticmethod
    def digest(counts, declared):
        key = (tuple(sorted(counts.items())), tuple(sorted(declared.items())))
        return hashlib.blake2b(repr(key).encode(), digest_size=16).digest()

    def grow(self, seed, fout):
        """Branch on half-edge demands and on *tight* slot demands (at most kmax exact satisfiers);
        slot demands with more satisfiers are left to the combine stage (their partners may lie
        in other clusters). With X > 0 a tight slot demand may also be declared residual (met
        inexactly in the member; at most 2X declared slots per cluster)."""
        empty = collections.Counter()
        a0, sv0, slots0, ex0, hh0, st0 = self.analyse({seed.id: 1}, empty)
        if not self.bounds_ok(a0, sv0, slots0, st0, 0):
            self.stats["dead_seed"] += 1
            return
        stack = [(collections.Counter({seed.id: 1}), empty)]
        K = self.kmax
        while stack:
            counts, declared = stack.pop()
            dig = self.digest(counts, declared)
            if dig in self.visited:
                self.stats["revisit"] += 1
                continue
            self.visited.add(dig)
            self.stats["states"] += 1
            a, sv, slots, ex, hh, st = self.analyse(counts, declared)
            ndecl = sum(declared.values())
            demands = []
            for t in range(NT):
                if a[t][t] % 2:
                    demands.append([(w, None) for w in self.by_half.get((t, t), ()) if w.id >= seed.id])
                for m in range(NT):
                    if m != t and a[t][m] > a[m][t]:
                        demands.append([(w, None) for w in self.by_half.get((m, t), ()) if w.id >= seed.id])
            nres = 0
            for j in range(3):
                for fv, cnt in sv[j].items():
                    if fv == (2, 2) or cnt <= 0:
                        continue
                    other = sv[j].get((4 - fv[0], 4 - fv[1]), 0)
                    if cnt <= other:
                        continue
                    comp = (4 - fv[0], 4 - fv[1])
                    opts = [(w, None) for w in self.by_slot.get((j, comp), ()) if w.id >= seed.id] if comp[0] >= 0 and comp[1] >= 0 else []
                    if len(opts) > K:
                        nres += cnt - other
                        continue
                    if ndecl < 2 * self.X:
                        opts.append((None, (j, fv)))
                    demands.append(opts)
            if not demands:
                self.record(counts, declared, sv, slots, ex, hh, st, nres, fout)
                continue
            cand = min(demands, key=len)
            if not cand:
                self.stats["dead_end"] += 1
                continue
            for (w, decl) in cand:
                self.stats["children"] += 1
                used2 = ndecl + (1 if decl else 0)
                if w is not None:
                    self.add(w, 1, a, sv, slots, ex, hh, st)
                if decl:
                    sv[decl[0]][decl[1]] -= 1
                ok = self.bounds_ok(a, sv, slots, st, (used2 + 1) // 2)
                if decl:
                    sv[decl[0]][decl[1]] += 1
                if w is not None:
                    self.add(w, -1, a, sv, slots, ex, hh, st)
                if not ok:
                    self.stats["pruned_child"] += 1
                    continue
                c2 = counts.copy()
                if w is not None:
                    c2[w.id] += 1
                d2 = declared
                if decl:
                    d2 = declared.copy()
                    d2[decl] += 1
                if self.digest(c2, d2) in self.visited:
                    self.stats["revisit"] += 1
                    continue
                stack.append((c2, d2))

    def record(self, counts, declared, sv, slots, ex, hh, st, nres, fout):
        self.stats["clusters"] += 1
        ndecl = sum(declared.values())
        self.stats[f"clusters_decl{ndecl}"] += 1
        rec = {"W": sorted(counts.items()), "deficit": st["deficit"], "cost8": st["cost8"], "slots": slots, "G": st["G"],
               "quads": st["quads"], "ex": ex, "hh": hh, "nres": nres + ndecl, "decl": ndecl,
               "sv22": [sv[j].get((2, 2), 0) for j in range(3)]}
        fout.write(json.dumps(rec) + "\n")


def stage_clusters(args):
    t0 = time.time()
    N = int(args[0])
    cmax = FR(next(a.split("=")[1] for a in args if a.startswith("--cmax=")))
    X = int(next(a.split("=")[1] for a in args if a.startswith("--excess=")))
    maxsize = int(next((a.split("=")[1] for a in args if a.startswith("--maxsize=")), N))
    out = next(a.split("=")[1] for a in args if a.startswith("--out="))
    logp = next((a.split("=")[1] for a in args if a.startswith("--log=")), None)
    log = open(logp, "a") if logp else sys.stderr
    closure = "--closure" in args
    cmax8 = int(cmax * 8)
    wins = window_table(cmax8, closure, log)
    kmax = int(next((a.split("=")[1] for a in args if a.startswith("--kmax=")), 3))
    g = Growth(wins, N, cmax8, X, maxsize, log, kmax=kmax)
    g.grow_all(out)
    summary = {"N": N, "cmax": str(cmax), "excess": X, "kmax": kmax, "maxsize": maxsize, "closure": closure, "windows": len(wins),
               "stats": dict(g.stats), "visited": len(g.visited), "wall_s": round(time.time() - t0, 1), "out": out}
    print(json.dumps(summary), file=log, flush=True)
    # window table for the combine stage
    json.dump({"summary": summary, "windows": [[list(w.x), [list(t) for t in w.M]] for w in wins]}, open(out + ".meta.json", "w"))


def sigma_of(c, wtab):
    """Signed slot imbalance vector of a cluster: for every block j and complementary pair
    {a, 4 - a} (a < 4 - a), count(a) - count(4 - a); (2,2) is self-complementary and omitted."""
    svc = [collections.Counter() for _ in range(3)]
    for wid, cnt in c["W"]:
        for (j, fv) in wtab[wid].svec:
            svc[j][fv] += cnt
    sig = {}
    for j in range(3):
        for fv, cnt in svc[j].items():
            if fv == (2, 2) or cnt == 0:
                continue
            b = (4 - fv[0], 4 - fv[1])
            if fv < b:
                sig[(j, fv)] = sig.get((j, fv), 0) + cnt
            else:
                sig[(j, b)] = sig.get((j, b), 0) - cnt
    return {k: v for k, v in sig.items() if v}, [svc[j].get((2, 2), 0) for j in range(3)]


def stage_combine(args):
    """Multisets of clusters forming a shape. Demand-driven at the cluster level: the union's slot
    imbalance ||sum sigma||_1 must end <= 2 excess (0 for tight shapes); a coordinate that exceeds
    the allowance must be reduced by a cluster of opposite sign there, so the search branches on the
    coordinate with the fewest such clusters; when nothing is forced (imbalance within the
    allowance) any cluster of id >= the seed may be added. Memoised by multiset digest."""
    t0 = time.time()
    N = int(args[0])
    inp = next(a.split("=")[1] for a in args if a.startswith("--clusters="))
    accept = next(a.split("=")[1] for a in args if a.startswith("--accept="))
    accept8 = {int(FR(v) * 8) for v in accept.split(",")}
    cmax8 = max(accept8)
    out = next(a.split("=")[1] for a in args if a.startswith("--out="))
    logp = next((a.split("=")[1] for a in args if a.startswith("--log=")), None)
    log = open(logp, "a") if logp else sys.stderr
    Dtot = 3 * N - 96
    excess_max = Dtot - min(accept8) // 8
    A = 2 * excess_max
    meta = json.load(open(inp + ".meta.json"))
    wtab = [W.Window(k, tuple(x), tuple(tuple(t) for t in M)) for k, (x, M) in enumerate(meta["windows"])]
    wid_of = {(w.x, w.M): w.id for w in wtab}
    cl = [json.loads(l) for l in open(inp) if l.strip()]
    cl = [c for c in cl if c["cost8"] <= cmax8 and c["decl"] <= A]
    for c in cl:
        c["sigma"], c["sv22"] = sigma_of(c, wtab)
        c["l1"] = sum(abs(v) for v in c["sigma"].values())
    # a cluster whose own imbalance cannot be repaired by the rest (deficit, cost) is useless:
    # the other clusters add at most 2 (Dtot - deficit) + 3 * 4 * (cmax - cost) slots ... keep simple
    print(f"{len(cl)} clusters loaded (cost <= {cmax8 / 8}), {time.time() - t0:.0f} s", file=log, flush=True)
    cl.sort(key=lambda c: (c["l1"] == 0, -c["deficit"], c["cost8"]))
    L = len(cl)
    for k, c in enumerate(cl):
        c["id"] = k
    fixers = collections.defaultdict(list)   # (key, sign): clusters with sigma[key] of that sign
    for c in cl:
        for key, v in c["sigma"].items():
            fixers[(key, 1 if v > 0 else -1)].append(c)
    stats = collections.Counter()
    seen = set()
    visited = set()
    fout = open(out, "w")
    MINTRIP8 = min((w.cost8 for w in wtab if w.nx == 3), default=2)

    def addable(c, tot, seedid):
        if c["id"] < seedid:
            return False
        nc = tot["cost8"] + c["cost8"]
        if nc > cmax8:
            return False
        nd = tot["deficit"] + c["deficit"]
        if (c["deficit"] >= 0 and nd > Dtot + 2) or nd < -2:
            return False
        if tot["G"] + c["G"] > N:
            return False
        if any(tot["slots"][j] + c["slots"][j] > 32 for j in range(3)):
            return False
        if tot["decl"] + c["decl"] > 2 * (Dtot - (nc + 7) // 8):
            return False
        return True

    def totals(chosen):
        tot = {"cost8": 0, "deficit": 0, "G": 0, "slots": [0, 0, 0], "decl": 0, "quads": 0}
        sig = collections.Counter()
        for k, m in chosen.items():
            c = cl[k]
            tot["cost8"] += m * c["cost8"]
            tot["deficit"] += m * c["deficit"]
            tot["G"] += m * c["G"]
            tot["decl"] += m * c["decl"]
            tot["quads"] += m * c["quads"]
            for j in range(3):
                tot["slots"][j] += m * c["slots"][j]
            for key, v in c["sigma"].items():
                sig[key] += m * v
        return tot, sig

    def complete(chosen, tot, sig):
        cost8 = tot["cost8"]
        if cost8 not in accept8 or tot["deficit"] != Dtot:
            return
        G = tot["G"]
        if (N - G) % 2:
            return
        u = (N - G) // 2
        g = 32 - 2 * u
        if tot["slots"] != [g, g, g]:
            stats["unbalanced"] += 1
            return
        stats["balanced"] += 1
        excess = Dtot - cost8 // 8
        if sum(abs(v) for v in sig.values()) > 2 * excess:
            stats["imbalance_fail"] += 1
            return
        ex = [sum(m * cl[k]["ex"][i] for k, m in chosen.items()) for i in range(6)]
        if any(v < 0 for v in ex):
            stats["rows_fail"] += 1
            return
        hh = [sum(m * cl[k]["hh"][i] for k, m in chosen.items()) for i in range(3)]
        if any(hh[i] + 2 * u < 16 for i in range(3)):
            stats["hubhub_fail"] += 1
            return
        if excess == 0:
            for j in range(3):
                if (sum(m * cl[k]["sv22"][j] for k, m in chosen.items()) + 2 * u) % 2:
                    stats["parity_fail"] += 1
                    return
        counts = collections.Counter()
        for k, m in chosen.items():
            for wid, cnt in cl[k]["W"]:
                counts[wid] += m * cnt
        items = sorted(counts.items())
        key = W.canon_W(items, wtab, wid_of)
        dig = hashlib.blake2b(repr(key).encode(), digest_size=16).digest()
        if dig in seen:
            stats["relabel_dup"] += 1
            return
        seen.add(dig)
        svl = []
        for wid, cnt in items:
            svl.extend(wtab[wid].svec * cnt)
        if not S.matching_ok([W.sv_gadget(svl)]):
            stats["sv_matching_fail"] += 1
            return
        stats["survive"] += 1
        stats[f"survive_cost_{cost8 // 8}"] += 1
        rec_ = {"u": u, "G": G, "cost": cost8 // 8, "phi": str(-FR(cost8, 8)), "clusters": sum(chosen.values()),
                "W": [[list(wtab[wid].x), [list(t) for t in wtab[wid].M], cnt] for wid, cnt in items]}
        fout.write(json.dumps(rec_) + "\n")
        fout.flush()

    for seed in cl:
        stack = [collections.Counter({seed["id"]: 1})]
        while stack:
            chosen = stack.pop()
            dig = hashlib.blake2b(repr(tuple(sorted(chosen.items()))).encode(), digest_size=16).digest()
            if dig in visited:
                continue
            visited.add(dig)
            stats["nodes"] += 1
            tot, sig = totals(chosen)
            l1 = sum(abs(v) for v in sig.values())
            if tot["deficit"] == Dtot:
                complete(chosen, tot, sig)
            # the remaining clusters can add at most this many slots (bounds the repairable imbalance)
            R = Dtot - tot["deficit"]
            C = cmax8 - tot["cost8"]
            if R < 0:
                R = 0
            maxslots = 2 * R + 3 * (C // MINTRIP8) + 4 * (C // 40 if C >= 40 else 0)
            if l1 > A + maxslots:
                stats["pruned_imbalance"] += 1
                continue
            if tot["deficit"] >= Dtot + 2 and l1 == 0:
                continue
            cand = None
            if l1 > A:
                # some coordinate must be reduced; the one over the allowance with the fewest fixers,
                # else the union of fixers of all nonzero coordinates
                over = [(key, v) for key, v in sig.items() if abs(v) > A]
                keys = over if over else [(key, v) for key, v in sig.items() if v]
                best = None
                for key, v in keys:
                    opts = [c for c in fixers.get((key, -1 if v > 0 else 1), ()) if addable(c, tot, seed["id"])]
                    if over:
                        if best is None or len(opts) < len(best):
                            best = opts
                    else:
                        best = (best or []) + opts
                cand = {c["id"]: c for c in (best or [])}.values()
            elif tot["deficit"] < Dtot or tot["cost8"] not in accept8 or True:
                # nothing forced: any cluster of id >= seed (deficit still to be filled, or a
                # deficit-0 cluster such as a triple pair)
                if tot["deficit"] < Dtot + 2:
                    cand = [c for c in cl if addable(c, tot, seed["id"])]
                else:
                    cand = []
            for c in cand:
                c2 = chosen.copy()
                c2[c["id"]] += 1
                stack.append(c2)
        stats["seeds"] += 1
        if stats["seeds"] % 20000 == 0:
            print(f"seed {stats['seeds']}/{L}: {json.dumps(dict(stats))}, {time.time() - t0:.0f} s", file=log, flush=True)
    fout.close()
    summary = {"N": N, "accept": accept, "clusters": L, "stats": dict(stats), "wall_s": round(time.time() - t0, 1), "out": out}
    print(json.dumps(summary), file=log, flush=True)


if __name__ == "__main__":
    {"clusters": stage_clusters, "combine": stage_combine}[sys.argv[1]](sys.argv[2:])
