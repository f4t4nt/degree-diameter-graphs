"""Window-multiset enumeration of shapes for three K4,4 blocks at n = 72 + N (TEMPLATE_CALCULUS.md 13.7).

A *window* is a connector type x (|x| <= 3, degree 4 - |x|) together with the multiset M of the
types of its connector neighbours; quads (|x| = 4) are windows with M empty. Every shape quantity
is a sum over windows: f = x + sum(M), phi = |x||f| - x.f - 3|x| - 3, deficit 3 - |x| (quads -1),
slots x, rows x_i (f_j - 2), and the slot vectors (block j, (f_i)_{i != j}) of the x_j slots in
block j. The potential Phi = sum over connector edges of E(x, y) <= 0 is written as a sum of
non-negative window costs cost(w) = -sum_{t in M} E(x, t) / 2 (each edge is paid half by each end),
so -Phi = sum_w cost(w) is a bounded knapsack: at n = 108 only windows of cost <= 12 can occur.

A gadget window is any window that is not a unit window (a triple whose neighbour is its
complement, x + t = (2,2,2)); gadget windows never touch unit members, so the gadget part of a
member is a multiset W of gadget windows with: deficit sum D = 3N - 96 (quads count -1), cost sum
= -Phi, |W| = G with u = (N - G) / 2 units, slots (32 - 2u) x 3, the six rows >= 0, the chaining
(Eulerian) condition a(t -> m) = a(m -> t) on half-edge counts (a(t -> t) even), the port lemma
(inequality form) as a slot matching, hub-to-hub coverage.

Enumeration: depth-first over the windows in a fixed order (types in a chosen order, quads last,
cost increasing inside a type) with multiplicities. Pruning (all necessary conditions):
  - cost + chaining bound <= Cmax: every unmatched half-edge (t -> m) will be paid for by a
    partner half-edge (m -> t) of cost -E(t, m) / 2;
  - cost + slot-complement bound <= Cmax: a slot (j, a) without a present complement (j, 4 - a)
    needs a future window providing one, paid at min cost(w') / |x(w')| per slot; the tolerance
    2 (D - cost) of most expensive unmatched slots is dropped (the total excess Phi + D bounds the
    number of inexact (port, block) pairs, each involving two slots);
  - deficit <= D + affordable quads; demanded types must fit the remaining deficit;
  - half-edge counts become final when the last window carrying them is passed: equality with the
    partner count if that is final too, otherwise an upper bound; a half-edge (t -> m) whose
    partner count (m -> t) is already final must not exceed it.
Every complete multiset is checked exactly: deficit, cost in the accepted set, slot balance, rows,
G parity, hub-to-hub, (2,2)-slot parity, slot-matching relaxation (one connector per slot).
Survivors are canonicalised under the S3 of block relabellings and streamed to a jsonl file.

Stage shapes: a window multiset is realised as all simple typed graphs on G vertices with the
prescribed neighbour-type multisets (depth-first, first-fresh symmetry breaking, pynauty
certificates); each realisation is split into connected components (gadgets), checked with the
port-lemma matching of shapes107 (inequality form) and the excess-bounded matching below (exact
complements with (b) and (a') at exact coordinates, at most Phi + D excess in total); shapes are
deduplicated under block relabelling and written in the format of shapes108.py.

Usage:
  windows108.py enum N --cmax=C [--accept=10,12] --out=runs/W.jsonl [--log=FILE] [--order=NAME]
  windows108.py shapes --in=runs/W.jsonl --out=runs/shapes.json [--log=FILE] [--oldexact] [--noexcess]
  windows108.py compare NEW.json OLD.json ...      (canonical shape keys, set comparison)
"""
import collections
import hashlib
import itertools
import json
import pathlib
import sys
import time
from fractions import Fraction as FR

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import components107 as C  # noqa: E402
import shapes107 as S  # noqa: E402
import shapes_fast as SF  # noqa: E402

NONISO = C.NONISO            # 20 types with |x| <= 3
QUADS = C.QUADS              # 15 types with |x| = 4
TIDX = {t: i for i, t in enumerate(NONISO)}
ETAB = C.ETAB                # E(x, y) as Fractions, all <= 0
PAIRS = S.PAIRS              # six ordered block pairs
PERMS = list(itertools.permutations(range(3)))
BIG = 10 ** 9


def is_unit_window(x, M):
    return sum(x) == 3 and len(M) == 1 and tuple(a + b for a, b in zip(x, M[0])) == (2, 2, 2)


class Window:
    __slots__ = ("id", "x", "M", "f", "cost8", "deficit", "phi", "ex", "hh", "svec", "halves", "nx", "quad", "t")

    def __init__(self, wid, x, M):
        self.id = wid
        self.x = x
        self.M = tuple(sorted(M))
        self.nx = sum(x)
        self.quad = self.nx == 4
        self.t = TIDX[x] if not self.quad else None
        f = list(x)
        for t in M:
            for j in range(3):
                f[j] += t[j]
        self.f = tuple(f)
        if self.quad:
            cost = FR(15 - C.self_w(x))
            self.deficit = -1
        else:
            cost = -sum(ETAB[(x, t)] for t in M) / 2
            self.deficit = 3 - self.nx
        assert cost >= 0 and (cost * 8).denominator == 1, (x, M, cost)
        self.cost8 = int(cost * 8)
        self.phi = self.nx * sum(f) - sum(a * b for a, b in zip(x, f)) - 3 * self.nx - 3
        self.ex = tuple(x[i] * (f[j] - 2) for (i, j) in PAIRS)
        self.hh = (x[0] * x[1], x[0] * x[2], x[1] * x[2])
        self.svec = []  # (block j, foreign vector) with multiplicity x_j
        for j in range(3):
            fv = tuple(f[i] for i in range(3) if i != j)
            for _ in range(x[j]):
                self.svec.append((j, fv))
        self.halves = tuple(sorted(collections.Counter(TIDX[t] for t in M).items()))


def all_windows(cmax8):
    wins = []
    for x in NONISO:
        d = 4 - sum(x)
        for M in itertools.combinations_with_replacement(NONISO, d):
            if is_unit_window(x, M):
                continue
            w = Window(len(wins), x, M)
            if w.cost8 <= cmax8:
                wins.append(w)
    for x in QUADS:
        w = Window(len(wins), x, ())
        if w.cost8 <= cmax8:
            wins.append(w)
    return wins


def type_order(name):
    trip = [t for t in NONISO if sum(t) == 3]
    dbl = [t for t in NONISO if sum(t) == 2]
    sgl = [t for t in NONISO if sum(t) == 1]
    rtr = [t for t in NONISO if sum(t) == 0]
    groups = {"t": trip, "d": dbl, "s": sgl, "r": rtr}
    out = []
    for ch in name:
        out += groups[ch]
    assert len(out) == 20, name
    return out


def canon_W(items, wins, wid_of):
    """Canonical key of a window multiset under block relabelling: min over S3 of the sorted list
    of (window id of the relabelled window, count)."""
    best = None
    for p in PERMS:
        key = []
        for wid, cnt in items:
            w = wins[wid]
            x = tuple(w.x[p[j]] for j in range(3))
            M = tuple(sorted(tuple(t[p[j]] for j in range(3)) for t in w.M))
            key.append((wid_of[(x, M)], cnt))
        key.sort()
        key = tuple(key)
        if best is None or key < best:
            best = key
    return best


def sv_gadget(svl):
    """Synthetic gadget with one connector per slot (the relaxation of shapes108.sv_gadget)."""
    types, f, slots = [], [], []
    for h, (j, fv) in enumerate(svl):
        t = [0, 0, 0]
        t[j] = 1
        types.append(tuple(t))
        full = list(fv)
        full.insert(j, 0)
        f.append(tuple(full))
        slots.append((h, j, tuple(4 - full[i] if i != j else 0 for i in range(3))))
    return {"types": types, "f": f, "slots_list": slots, "edges": []}


# ----------------------------------------------------------------------------------------------
# stage enum
# ----------------------------------------------------------------------------------------------

def stage_enum(args):
    t0 = time.time()
    N = int(args[0])
    cmax = FR(next(a.split("=")[1] for a in args if a.startswith("--cmax=")))
    accept = next((a.split("=")[1] for a in args if a.startswith("--accept=")), None)
    accept8 = {int(FR(v) * 8) for v in accept.split(",")} if accept else None
    out = next((a.split("=")[1] for a in args if a.startswith("--out=")), None)
    logp = next((a.split("=")[1] for a in args if a.startswith("--log=")), None)
    order = next((a.split("=")[1] for a in args if a.startswith("--order=")), "srdt")
    nodelimit = int(next((a.split("=")[1] for a in args if a.startswith("--nodelimit=")), 0))
    log = open(logp, "a") if logp else sys.stderr
    cmax8 = int(cmax * 8)
    Dtot = 3 * N - 96
    assert Dtot >= 0
    wins = all_windows(cmax8)
    wid_of = {(w.x, w.M): w.id for w in wins}
    torder = {t: k for k, t in enumerate(type_order(order))}
    seq = sorted([w for w in wins if not w.quad], key=lambda w: (torder[w.x], w.cost8, w.M)) + \
          sorted([w for w in wins if w.quad], key=lambda w: w.cost8)
    L = len(seq)
    print(f"N={N} Dtot={Dtot} cmax={cmax} windows={len(wins)} (quads {sum(w.quad for w in wins)}) order={order}", file=log, flush=True)
    print(f"windows by |x|: {dict(sorted(collections.Counter(w.nx for w in wins).items()))}", file=log, flush=True)
    NT = len(NONISO)
    DEG = [4 - sum(t) for t in NONISO]
    DEF = [3 - sum(t) for t in NONISO]
    HC = [[int(-ETAB[(NONISO[t], NONISO[m])] / 2 * 8) for m in range(NT)] for t in range(NT)]
    # type blocks in seq: block_start[k] = first position of the next type block after position k
    next_block = [L] * (L + 1)
    k = L - 1
    cur = None
    while k >= 0:
        key = seq[k].t if not seq[k].quad else "q"
        if cur is None or key != cur:
            cur = key
            nb = k + 1
        next_block[k] = nb
        k -= 1
    # finalisation of half-edge counts
    last_pos = {}
    for k, w in enumerate(seq):
        if w.quad:
            continue
        for m, _ in w.halves:
            last_pos[(w.t, m)] = k
    final_at = collections.defaultdict(list)
    for (t, m), k in last_pos.items():
        final_at[k + 1].append((t, m))
    # min cost per slot of a window providing slot vector (j, b) (eighths, fractional allowed)
    minsc = {}
    for w in wins:
        if w.nx == 0:
            continue
        c = w.cost8 / w.nx
        for (j, fv) in w.svec:
            if (j, fv) not in minsc or c < minsc[(j, fv)]:
                minsc[(j, fv)] = c

    def partner_cost(j, a):
        b = (4 - a[0], 4 - a[1])
        if b[0] < 0 or b[1] < 0:
            return BIG
        return minsc.get((j, b), BIG)

    MINTRIP8 = min(w.cost8 for w in wins if w.nx == 3)
    quad_costs = sorted(w.cost8 for w in wins if w.quad)
    cheapest_quad8 = quad_costs[0] if quad_costs else BIG

    # state
    a = [[0] * NT for _ in range(NT)]
    sv = [collections.Counter() for _ in range(3)]
    slots = [0, 0, 0]
    ex = [0] * 6
    hh = [0, 0, 0]
    st = {"deficit": 0, "cost8": 0, "G": 0, "quads": 0, "chain8": 0}
    stats = collections.Counter()
    chosen = []
    seen = set()
    fout = open(out, "w") if out else None
    results = []

    def add(w, sign):
        t = w.t
        if t is not None:
            for m, cnt in w.halves:
                for _ in range(cnt):
                    if sign > 0:
                        if t == m:
                            st["chain8"] += HC[t][t] if a[t][t] % 2 == 0 else -HC[t][t]
                            a[t][t] += 1
                        else:
                            st["chain8"] += HC[t][m] if a[t][m] >= a[m][t] else -HC[t][m]
                            a[t][m] += 1
                    else:
                        if t == m:
                            a[t][t] -= 1
                            st["chain8"] -= HC[t][t] if a[t][t] % 2 == 0 else -HC[t][t]
                        else:
                            a[t][m] -= 1
                            st["chain8"] -= HC[t][m] if a[t][m] >= a[m][t] else -HC[t][m]
        for (j, fv) in w.svec:
            sv[j][fv] += sign
        for j in range(3):
            slots[j] += sign * w.x[j]
            hh[j] += sign * w.hh[j]
        for k in range(6):
            ex[k] += sign * w.ex[k]
        st["deficit"] += sign * w.deficit
        st["cost8"] += sign * w.cost8
        st["G"] += sign
        if w.quad:
            st["quads"] += sign

    def slot_bound_ok():
        cost8 = st["cost8"]
        budget = cmax8 - cost8
        tol = 2 * (Dtot - (cost8 + 7) // 8)
        if tol < 0:
            return False
        costs = []
        for j in range(3):
            c = sv[j]
            for fv, cnt in c.items():
                if cnt <= 0 or fv == (2, 2):
                    continue
                other = c.get((4 - fv[0], 4 - fv[1]), 0)
                if cnt > other:
                    costs.extend([partner_cost(j, fv)] * (cnt - other))
        if len(costs) <= tol:
            return True
        costs.sort(reverse=True)
        return sum(costs[tol:]) <= budget + 1e-9

    def demand_ok():
        qmax = (cmax8 - st["cost8"]) // cheapest_quad8 if cheapest_quad8 < BIG else 0
        R = Dtot - st["deficit"] + qmax
        need_def = 0
        for m in range(NT):
            if DEF[m] == 0:
                continue
            D = a[m][m] % 2
            am = a[m]
            for t in range(NT):
                if t != m:
                    d = a[t][m] - am[t]
                    if d > 0:
                        D += d
            if D:
                need_def += DEF[m] * (-(-D // DEG[m]))
                if need_def > R:
                    return False
        return True

    def finalise_ok(k):
        for (t, m) in final_at.get(k, ()):
            if t == m:
                if a[t][t] % 2:
                    return False
            elif last_pos.get((m, t), -1) < k:
                if a[t][m] != a[m][t]:
                    return False
            elif a[m][t] > a[t][m]:
                return False
        return True

    def partner_final_ok(w, k):
        """Half-edges (t -> m) of w whose partner count (m -> t) is final must not exceed it."""
        t = w.t
        for m, _ in w.halves:
            if last_pos.get((m, t), -1) < k and a[t][m] > a[m][t]:
                return False
        return True

    def complete():
        cost8 = st["cost8"]
        if accept8 is not None and cost8 not in accept8:
            return
        if cost8 % 8:
            return
        G = st["G"]
        if (N - G) % 2:
            return
        u = (N - G) // 2
        g = 32 - 2 * u
        if slots != [g, g, g]:
            return
        for t in range(NT):
            if a[t][t] % 2:
                return
            for m in range(t + 1, NT):
                if a[t][m] != a[m][t]:
                    return
        stats["complete_chain"] += 1
        if any(v < 0 for v in ex):
            stats["rows_fail"] += 1
            return
        if any(hh[k] + 2 * u < 16 for k in range(3)):
            stats["hubhub_fail"] += 1
            return
        for j in range(3):
            if (sv[j].get((2, 2), 0) + 2 * u) % 2:
                stats["parity_fail"] += 1
                return
        items = sorted(collections.Counter(chosen).items())
        key = canon_W(items, wins, wid_of)
        dig = hashlib.blake2b(repr(key).encode(), digest_size=16).digest()
        if dig in seen:
            stats["relabel_dup"] += 1
            return
        seen.add(dig)
        svl = []
        for wid, cnt in items:
            svl.extend(wins[wid].svec * cnt)
        if not S.matching_ok([sv_gadget(svl)]):
            stats["sv_matching_fail"] += 1
            return
        stats["survive"] += 1
        stats[f"survive_cost_{cost8 // 8}"] += 1
        rec_ = {"u": u, "G": G, "cost": cost8 // 8, "phi": str(-FR(cost8, 8)),
                "W": [[list(wins[wid].x), [list(t) for t in wins[wid].M], cnt] for wid, cnt in items]}
        if fout:
            fout.write(json.dumps(rec_) + "\n")
        else:
            results.append(rec_)

    def rec(start):
        stats["nodes"] += 1
        if nodelimit and stats["nodes"] > nodelimit:
            raise RuntimeError("node limit")
        if st["deficit"] == Dtot and chosen:
            complete()
        if stats["nodes"] % 1000000 == 0:
            print(f"nodes {stats['nodes']}: {json.dumps(dict(stats))}, {time.time() - t0:.0f} s", file=log, flush=True)
        k = start
        while k < L:
            if not finalise_ok(k):
                break
            w = seq[k]
            cost8 = st["cost8"]
            nc = cost8 + w.cost8
            if nc > cmax8:
                k = next_block[k]  # windows of this type are sorted by cost
                continue
            qmax = (cmax8 - nc) // cheapest_quad8 if cheapest_quad8 < BIG else 0
            nd = st["deficit"] + w.deficit
            if not w.quad and nd > Dtot + qmax:
                if w.deficit > 0:
                    k = next_block[k]  # every window of this type has the same deficit
                    continue
                k += 1
                continue
            if w.quad and nd < Dtot:
                k += 1
                continue
            if slots[0] + w.x[0] > 32 or slots[1] + w.x[1] > 32 or slots[2] + w.x[2] > 32:
                k += 1
                continue
            add(w, +1)
            ok = st["chain8"] + st["cost8"] <= cmax8
            if ok and not w.quad:
                ok = partner_final_ok(w, k)
            if ok:
                ok = demand_ok()
            if ok:
                ok = slot_bound_ok()
            if ok:
                R = Dtot - st["deficit"]
                qm = (cmax8 - st["cost8"]) // cheapest_quad8 if cheapest_quad8 < BIG else 0
                need_fill = 3 * max(slots) - sum(slots)
                if need_fill > 2 * max(R + qm, 0) + 3 * ((cmax8 - st["cost8"]) // MINTRIP8) + 4 * qm:
                    ok = False
            if ok:
                chosen.append(w.id)
                rec(k)
                chosen.pop()
            else:
                stats["pruned"] += 1
            add(w, -1)
            k += 1

    try:
        rec(0)
    except RuntimeError as e:
        stats["aborted"] = str(e)
    if fout:
        fout.close()
    summary = {"N": N, "cmax": str(cmax), "accept": accept, "order": order, "windows": len(wins), "seq": L,
               "stats": dict(stats), "wall_s": round(time.time() - t0, 1), "out": out}
    print(json.dumps(summary), file=log, flush=True)
    if not out:
        json.dump({"summary": summary, "W": results}, sys.stdout)
    return summary


# ----------------------------------------------------------------------------------------------
# stage shapes: realisation, components, matching, dedup
# ----------------------------------------------------------------------------------------------

def realisations(wlist):
    """All simple typed graphs on the connectors wlist (one Window per connector) in which every
    connector's neighbour-type multiset is its window, up to isomorphism (typed certificates).
    Depth-first: the first connector with open needs is completed type by type, partners in
    increasing index order; among fresh (isolated so far) connectors with the same window only the
    first is tried (they are interchangeable)."""
    G = len(wlist)
    need = [collections.Counter(w.M) for w in wlist]
    adj = [set() for _ in range(G)]
    edges = []
    found = {}
    types = [w.x for w in wlist]
    last = {}

    def rec():
        v = next((i for i in range(G) if need[i]), None)
        if v is None:
            key = C.closed_cert(types, set(edges))
            if key not in found:
                found[key] = sorted(edges)
            return
        tau = min(need[v])
        tv = types[v]
        lo = last.get((v, tau), -1) + 1
        seen_fresh = set()
        for w in range(lo, G):
            if w == v or w in adj[v] or types[w] != tau or need[w].get(tv, 0) == 0:
                continue
            if not adj[w]:
                if wlist[w].id in seen_fresh:
                    continue
                seen_fresh.add(wlist[w].id)
            adj[v].add(w)
            adj[w].add(v)
            edges.append((min(v, w), max(v, w)))
            need[v][tau] -= 1
            if need[v][tau] == 0:
                del need[v][tau]
            need[w][tv] -= 1
            if need[w][tv] == 0:
                del need[w][tv]
            old = last.get((v, tau))
            last[(v, tau)] = w
            rec()
            if old is None:
                del last[(v, tau)]
            else:
                last[(v, tau)] = old
            need[w][tv] += 1
            need[v][tau] += 1
            edges.pop()
            adj[v].discard(w)
            adj[w].discard(v)

    rec()
    return [(types, e) for e in found.values()]


def components_of(types, edges):
    n = len(types)
    nbr = {i: set() for i in range(n)}
    for a_, b_ in edges:
        nbr[a_].add(b_)
        nbr[b_].add(a_)
    seen = set()
    comps = []
    for s in range(n):
        if s in seen:
            continue
        comp = []
        stack = [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in nbr[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        comp.sort()
        idx = {v: k for k, v in enumerate(comp)}
        ctypes = [list(types[v]) for v in comp]
        cedges = sorted([idx[a_], idx[b_]] for a_, b_ in edges if a_ in idx and b_ in idx)
        pot = sum(ETAB[(types[a_], types[b_])] for a_, b_ in edges if a_ in idx and b_ in idx) if cedges else C.iso_pot(tuple(ctypes[0]))
        c = {"types": ctypes, "edges": cedges, "n": len(comp), "slots": [sum(t[j] for t in ctypes) for j in range(3)], "potential": str(pot)}
        comps.append(S.analyse(c))
    return comps


def matching_excess(gadgets, u, excess, old_a=False):
    """Port lemma as a matching of gadget slots with bounded total excess.

    A slot of connector h in block j has foreign vector a = (f_i(h))_{i != j}. Slots with a = (2,2)
    are free (unit partners). Every other slot is paired with a gadget slot of the same block on
    another connector with f-sums >= 4 (excess = sum of the surpluses), or, when a >= (2,2), with a
    unit slot (excess |a| - 4, at most 2u such pairings per block). At every exact coordinate i
    (f-sum exactly 4): (b) a common connector neighbour of the two holders holds no block-i port;
    (a') adjacent holders both hold no block-i port (old_a: adjacent holders are excluded, the
    107 form, sound while component potentials stay above -10). The total excess must not
    exceed Phi + D, which is the exact total excess of the shape. Exhaustive search."""
    allslots = []
    for gi, g in enumerate(gadgets):
        nbr = {i: set() for i in range(len(g["types"]))}
        for a_, b_ in g["edges"]:
            nbr[a_].add(b_)
            nbr[b_].add(a_)
        types = [tuple(t) for t in g["types"]]
        for (h, j, need) in g["slots_list"]:
            allslots.append(((gi, h), j, need, tuple(g["f"][h]), nbr[h], types))
    weak = [k for k, s in enumerate(allslots) if any(v != 2 for i, v in enumerate(s[2]) if i != s[1])]
    if not weak:
        return True
    n = len(allslots)

    def pair_cost(k, l):
        (c1, j1, need1, f1, n1, types), (c2, j2, need2, f2, n2, _) = allslots[k], allslots[l]
        if c1 == c2 or j1 != j2:
            return None
        cost = 0
        for i in range(3):
            if i == j1:
                continue
            s_ = f1[i] + f2[i]
            if s_ < 4:
                return None
            if s_ > 4:
                cost += s_ - 4
                continue
            if c1[0] == c2[0]:
                if c2[1] in n1:
                    if old_a or types[c1[1]][i] or types[c2[1]][i]:
                        return None
                for c in n1 & n2:
                    if types[c][i]:
                        return None
        return cost

    cands = {}
    for k in weak:
        opts = []
        for l in range(n):
            if l == k:
                continue
            c = pair_cost(k, l)
            if c is not None and c <= excess:
                opts.append((l, c))
        s = allslots[k]
        fv = [s[3][i] for i in range(3) if i != s[1]]
        if min(fv) >= 2:
            uc = sum(fv) - 4
            if uc <= excess:
                opts.append((-1, uc))  # unit partner
        if not opts:
            return False
        cands[k] = sorted(opts, key=lambda o: o[1])
    order = sorted(weak, key=lambda k: len(cands[k]))
    used = set()
    unit_used = [0, 0, 0]

    def rec(idx, tol):
        if idx == len(order):
            return True
        k = order[idx]
        if k in used:
            return rec(idx + 1, tol)
        j = allslots[k][1]
        for l, c in cands[k]:
            if c > tol:
                break
            if l == -1:
                if unit_used[j] >= 2 * u:
                    continue
                unit_used[j] += 1
                used.add(k)
                if rec(idx + 1, tol - c):
                    return True
                used.discard(k)
                unit_used[j] -= 1
                continue
            if l in used:
                continue
            used.update((k, l))
            if rec(idx + 1, tol - c):
                return True
            used.difference_update((k, l))
        return False

    return rec(0, excess)


def stage_shapes(args):
    t0 = time.time()
    inp = next(a.split("=")[1] for a in args if a.startswith("--in="))
    out = next(a.split("=")[1] for a in args if a.startswith("--out="))
    logp = next((a.split("=")[1] for a in args if a.startswith("--log=")), None)
    log = open(logp, "a") if logp else sys.stderr
    old_a = "--oldexact" in args
    noexcess = "--noexcess" in args   # apply the excess matching only to tight shapes (old pipelines)
    Dtot = None
    N = None
    stats = collections.Counter()
    shapes = []
    seen = set()
    wcache = {}
    nrec = 0
    for line in open(inp):
        if not line.strip():
            continue
        r = json.loads(line)
        nrec += 1
        u = r["u"]
        G = r["G"]
        if N is None:
            N = 2 * u + G
            Dtot = 3 * N - 96
        phi = FR(r["phi"])
        excess = int(phi) + Dtot
        wlist = []
        for x, M, cnt in r["W"]:
            key = (tuple(x), tuple(tuple(t) for t in M))
            if key not in wcache:
                wcache[key] = Window(len(wcache), key[0], key[1])
            wlist.extend([wcache[key]] * cnt)
        reals = realisations(wlist)
        stats["realisations"] += len(reals)
        if not reals:
            stats["unrealisable"] += 1
        for types, edges in reals:
            gadgets = components_of(types, edges)
            if any(S.is_unit(g) for g in gadgets):
                stats["unit_component"] += 1
                continue
            stats["candidates"] += 1
            if not S.matching_ok(gadgets):
                stats["ineq_matching_fail"] += 1
                continue
            if excess == 0 or not noexcess:
                if not matching_excess(gadgets, u, excess, old_a=old_a):
                    stats["excess_matching_fail"] += 1
                    stats[f"excess_matching_fail_phi{int(phi)}"] += 1
                    continue
            k = (u, SF.shape_key(gadgets))
            if k in seen:
                stats["relabel_dup"] += 1
                continue
            seen.add(k)
            stats["survive"] += 1
            shapes.append({"units": u, "potential": str(phi), "budget": 768 + 2 * excess,
                           "gadgets": [{"types": g["types"], "edges": g["edges"], "potential": g["potential"], "f": g["f"]} for g in gadgets]})
        if nrec % 2000 == 0:
            print(f"{nrec} multisets: {json.dumps(dict(stats))}, {time.time() - t0:.0f} s", file=log, flush=True)
    shapes.sort(key=lambda s_: (-FR(s_["potential"]), len(s_["gadgets"]), sum(len(g["types"]) for g in s_["gadgets"])))
    by = collections.Counter((s_["budget"], s_["units"]) for s_ in shapes)
    summary = {"N": N, "n": 72 + (N or 0), "multisets": nrec, "stats": dict(stats), "wall_s": round(time.time() - t0, 1),
               "by_budget_units": {f"{b}/{u}": c for (b, u), c in sorted(by.items())}}
    print(json.dumps(summary), file=log, flush=True)
    json.dump({"N": N, "n": 72 + (N or 0), "note": "windows108.py shapes", "summary": summary, "shapes": shapes}, open(out, "w"), indent=0)


def stage_compare(args):
    files = [a for a in args if not a.startswith("--")]
    sets = []
    for p in files:
        d = json.load(open(p))
        keys = {}
        for s_ in d["shapes"]:
            keys[(s_["units"], SF.shape_key(s_["gadgets"]))] = s_
        sets.append((p, keys))
    p0, k0 = sets[0]
    allold = {}
    for p, k in sets[1:]:
        allold.update(k)
    only_new = [k0[k] for k in k0.keys() - allold.keys()]
    only_old = [allold[k] for k in allold.keys() - k0.keys()]
    out = {"new": p0, "new_count": len(k0), "old_files": files[1:], "old_count": len(allold),
           "common": len(k0.keys() & allold.keys()), "n_only_new": len(only_new), "n_only_old": len(only_old),
           "only_new": only_new[:5], "only_old": only_old[:5]}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    {"enum": stage_enum, "shapes": stage_shapes, "compare": stage_compare}[sys.argv[1]](sys.argv[2:])
