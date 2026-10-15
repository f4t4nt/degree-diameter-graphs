"""Shape enumeration for three K4,4 blocks at n = 72 + N, sub-family of total potential >= F, in
two resumable stages so that every process stays under 60 minutes and 2 GB.

Mathematics as in shapes_fast.py (components107 enumeration, fine signatures,
coarse recursion with range tables, fine ex-vector search, port-lemma matching in inequality
form, block-relabelling dedup). Two changes make the matching phase tractable for the -8 family
at 108 (3.4M gadgets):

  1. Matching classes. The inequality matching (shapes107.matching_ok, exact=False) reads only,
     per connector h of a gadget, its profile x(h) (which blocks it has slots in) and its carried
     vector f(h). Gadgets with the same multiset MS = {(x(h), f(h))} are therefore
     indistinguishable for it, and MS also determines deficit, potential (T(c) = |x||F| - x.F),
     slots, n and the ex-vector, so MS refines the fine signature. Coarser still, the slot-vector
     multiset SV = {(block, f(h)) per slot} gives a relaxation of the matching (slots of one
     connector may then pair with each other), a necessary condition shared by every MS class
     with that SV. Stage "match" runs the relaxed matching once per multiset of SV classes
     (cached), the exact one once per multiset of MS classes inside the survivors, and expands
     only those to gadget multisets (which all pass), then deduplicates under block relabelling.
  2. Checkpointing and memory. Stage "enum" streams every gadget to a binary file (signature id,
     MS key, packed types and edges) and keeps in memory only the class index (signature -> MS
     key -> count) besides the enumerator's digest sets; it writes the index and the fine shapes
     to a pickle. Stage "match" runs the MS-level matching for a slice of the fine shapes first
     (no gadgets needed), then scans the gadget file once for the surviving classes only, expands
     them and deduplicates under block relabelling; it can be run in several processes.

The exact (tight) matching is not used here: the sub-family F > 96 - 3N contains no tight shape.
Tight shapes (F = 96 - 3N) must not be run through this script (see recheck_tight.py for the
status of the equality matching).

Usage:
  shapes108.py enum N --family=F --store=runs/store.pkl [--log=FILE]
  shapes108.py match --store=runs/store.pkl --out=runs/shapes.json [--slice=A:B] [--log=FILE]
  shapes108.py merge --out=runs/shapes.json PART.json ...   (concatenate slices, dedup, sort)
"""
import collections
import itertools
import json
import pickle
import struct
import sys
import time
from fractions import Fraction as FR

import pathlib

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import components107 as C  # noqa: E402
import shapes107 as S  # noqa: E402
import shapes_fast as SF  # noqa: E402

PAIRS = S.PAIRS


def ms_bytes(types, edges):
    """MS class key: sorted (x(h), f(h)) pairs, 6 bytes per connector."""
    n = len(types)
    nbr = {i: set() for i in range(n)}
    for a, b in edges:
        nbr[a].add(b)
        nbr[b].add(a)
    f = [tuple(types[i][j] + sum(types[k][j] for k in nbr[i]) for j in range(3)) for i in range(n)]
    return b"".join(bytes(t) + bytes(fv) for t, fv in sorted(zip(types, f)))


def ms_gadget(ms):
    """Synthetic gadget record for shapes107.matching_ok(exact=False) from an MS class key."""
    types = [tuple(ms[6 * i:6 * i + 3]) for i in range(len(ms) // 6)]
    f = [tuple(ms[6 * i + 3:6 * i + 6]) for i in range(len(ms) // 6)]
    slots = []
    for h in range(len(types)):
        for j in range(3):
            for _ in range(types[h][j]):
                need = tuple(4 - f[h][i] if i != j else 0 for i in range(3))
                slots.append((h, j, need))
    return {"types": types, "f": f, "slots_list": slots, "edges": []}


def sv_key(ms):
    """Slot-vector multiset of an MS class (connector identity dropped): sorted (block, f) per slot."""
    out = []
    for i in range(len(ms) // 6):
        t, fv = ms[6 * i:6 * i + 3], ms[6 * i + 3:6 * i + 6]
        for j in range(3):
            out += [(j, fv)] * t[j]
    return tuple(sorted(out))


def sv_gadget(sv):
    """Synthetic gadget with one connector per slot (the relaxation of the matching that ignores
    which slots belong to the same connector): a necessary condition for every MS class with
    this slot-vector multiset."""
    types, f, slots = [], [], []
    for h, (j, fv) in enumerate(sv):
        t = [0, 0, 0]
        t[j] = 1
        types.append(tuple(t))
        f.append(tuple(fv))
        slots.append((h, j, tuple(4 - fv[i] if i != j else 0 for i in range(3))))
    return {"types": types, "f": f, "slots_list": slots, "edges": []}


def stage_enum(args):
    t0 = time.time()
    N = int(args[0])
    family = int(next(a.split("=")[1] for a in args if a.startswith("--family=")))
    storep = next(a.split("=")[1] for a in args if a.startswith("--store="))
    logp = next((a.split("=")[1] for a in args if a.startswith("--log=")), None)
    log = open(logp, "a") if logp else sys.stderr
    bound4 = 4 * family
    Dtot = 3 * N - 96
    assert family > 96 - 3 * N, "tight family: use the exact matching pipeline"
    store = collections.defaultdict(collections.Counter)  # sig -> ms key (bytes) -> count
    sig_id = {}
    stats = collections.Counter()
    gadf = open(storep + ".gad", "wb")

    def sink(types, edges, pot):
        if S.is_unit({"n": len(types), "types": types}):
            stats["units"] += 1
            return
        sig = SF.analyse_raw(types, edges)
        if sig[1] < bound4:
            stats["below_bound"] += 1
            return
        ms = ms_bytes(types, edges)
        store[sig][ms] += 1
        if sig not in sig_id:
            sig_id[sig] = len(sig_id)
        tb, eb = SF.pack(types, edges)
        gadf.write(struct.pack("<IBBB", sig_id[sig], len(ms), len(tb), len(eb)) + ms + tb + eb)
        stats["gadgets"] += 1

    summary = C.enumerate_components(N, log=log, sink=sink, family_bound=family)
    gadf.close()
    print(json.dumps(summary), file=log, flush=True)
    n_ms = sum(len(v) for v in store.values())
    print(f"{stats['gadgets']} gadgets, {len(store)} fine signatures, {n_ms} matching classes, {time.time() - t0:.0f} s", file=log, flush=True)
    coarse = collections.defaultdict(list)
    for sig in store:
        coarse[sig[:4]].append(sig)
    keys = sorted(coarse, key=lambda k: (k[0], -k[1]))
    maxdef = summary["maxdef"]
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
        for p in range(P + 1):
            k = (R, p)
            if k in lo3 and lo3[k][0] != INF and all(lo3[k][j] <= need[j] <= hi3[k][j] for j in range(3)):
                return True
        return False

    last_quad = max([i for i, k in enumerate(keys) if k[0] < 0], default=-1)
    coarse_shapes = []
    target = (-Dtot // 3,) * 3

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
    fine_shapes = []
    for u, pot4, chosen in coarse_shapes:
        cnt = collections.Counter(chosen)
        slots = list(cnt.items())
        options = [[(sig[4], sig) for sig in coarse[key]] for key, m in slots]
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
    used = set(s for _, _, sigs in fine_shapes for s in sigs)
    index = {sig: dict(store[sig]) for sig in used}
    n_ms_used = sum(len(v) for v in index.values())
    n_g_used = sum(sum(v.values()) for v in index.values())
    print(f"{len(fine_shapes)} fine shapes use {len(used)} signatures, {n_ms_used} matching classes, {n_g_used} gadgets; "
          f"{time.time() - t0:.0f} s", file=log, flush=True)
    with open(storep, "wb") as f:
        pickle.dump({"N": N, "family": family, "bound4": bound4, "enumeration": summary, "stats": dict(stats),
                     "fine_signatures": len(store), "coarse_keys": len(keys), "matching_classes": n_ms,
                     "fine_shapes": fine_shapes, "index": index, "sig_id": {sig: sig_id[sig] for sig in used},
                     "wall_s": round(time.time() - t0, 1)}, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"stored {storep} (+ .gad), {time.time() - t0:.0f} s", file=log, flush=True)


def stage_match(args):
    t0 = time.time()
    storep = next(a.split("=")[1] for a in args if a.startswith("--store="))
    out = next(a.split("=")[1] for a in args if a.startswith("--out="))
    logp = next((a.split("=")[1] for a in args if a.startswith("--log=")), None)
    log = open(logp, "a") if logp else sys.stderr
    sl = next((a.split("=")[1] for a in args if a.startswith("--slice=")), None)
    data = pickle.load(open(storep, "rb"))
    N, bound4 = data["N"], data["bound4"]
    fine_shapes = data["fine_shapes"]
    index, sig_id = data["index"], data["sig_id"]
    A, B = (int(x) for x in sl.split(":")) if sl else (0, len(fine_shapes))
    B = min(B, len(fine_shapes))
    print(f"loaded {storep}: {len(fine_shapes)} fine shapes, slice {A}:{B}, {time.time() - t0:.0f} s", file=log, flush=True)
    stats = collections.Counter()
    ms_cache, sv_cache, sv_result = {}, {}, {}
    sv_of_sig = {}
    for s_ in index:
        d = collections.defaultdict(list)
        for ms in sorted(index[s_]):
            d[sv_key(ms)].append(ms)
        sv_of_sig[s_] = d
    print(f"{sum(len(d) for d in sv_of_sig.values())} slot-vector classes for {sum(len(v) for v in index.values())} matching classes",
          file=log, flush=True)
    survivors = []  # (fine idx, Counter of (sig, ms))
    needed = set()
    next_report = A + 2000
    for idx in range(A, B):
        u, pot4, sigs = fine_shapes[idx]
        cnt = collections.Counter(sigs)
        sig_list = list(cnt.items())
        # level 1: slot-vector classes (relaxed matching, cached by the multiset of sv keys)
        sv_pools = [list(itertools.combinations_with_replacement(sorted(sv_of_sig[s_]), m)) for s_, m in sig_list]
        for sv_combo in itertools.product(*sv_pools):
            svs = tuple(sorted(sv for part in sv_combo for sv in part))
            stats["sv_multisets"] += 1
            r = sv_result.get(svs)
            if r is None:
                gs = []
                for sv in svs:
                    if sv not in sv_cache:
                        sv_cache[sv] = sv_gadget(sv)
                    gs.append(sv_cache[sv])
                r = S.matching_ok(gs)
                sv_result[svs] = r
            if not r:
                stats["sv_matching_fail"] += 1
                continue
            stats["sv_matching_ok"] += 1
            # level 2: MS classes inside the chosen sv classes
            ms_pools = []
            for (s_, m), part in zip(sig_list, sv_combo):
                for sv, k in collections.Counter(part).items():
                    ms_pools.append([(s_, c) for c in itertools.combinations_with_replacement(sv_of_sig[s_][sv], k)])
            for ms_combo in itertools.product(*ms_pools):
                pairs = [(s_, ms) for s_, part in ms_combo for ms in part]
                stats["ms_multisets"] += 1
                gs = []
                for _, ms in pairs:
                    if ms not in ms_cache:
                        ms_cache[ms] = ms_gadget(ms)
                    gs.append(ms_cache[ms])
                if not S.matching_ok(gs):
                    stats["ms_matching_fail"] += 1
                    continue
                stats["ms_matching_ok"] += 1
                ms_cnt = collections.Counter(pairs)
                survivors.append((idx, ms_cnt))
                needed.update(ms_cnt)
        if len(sv_result) > 200000:
            sv_result.clear()
        if idx + 1 >= next_report:
            next_report += 2000
            print(f"fine shape {idx + 1}/{B}: {json.dumps(dict(stats))}, {time.time() - t0:.0f} s", file=log, flush=True)
    print(f"phase 1 done: {json.dumps(dict(stats))}, {len(survivors)} surviving class multisets, {len(needed)} classes needed, "
          f"{time.time() - t0:.0f} s", file=log, flush=True)
    # phase 2: scan the gadget file for the needed classes
    id_sig = {v: k for k, v in sig_id.items()}
    gadgets = collections.defaultdict(list)
    with open(storep + ".gad", "rb") as f:
        while True:
            head = f.read(7)
            if not head:
                break
            sid, lms, ltb, leb = struct.unpack("<IBBB", head)
            body = f.read(lms + ltb + leb)
            sig = id_sig.get(sid)
            if sig is None:
                continue
            ms = body[:lms]
            if (sig, ms) not in needed:
                continue
            gadgets[(sig, ms)].append(SF.gadget_record(body[lms:lms + ltb], body[lms + ltb:]))
    print(f"phase 2 done: {sum(len(v) for v in gadgets.values())} gadgets loaded, {time.time() - t0:.0f} s", file=log, flush=True)
    shapes = []
    seen_keys = set()
    for idx, ms_cnt in survivors:
        u, pot4, sigs = fine_shapes[idx]
        pools = [list(itertools.combinations_with_replacement(gadgets[key], m)) for key, m in ms_cnt.items()]
        for combo in itertools.product(*pools):
            gl = [g for part in combo for g in part]
            stats["gadget_multisets"] += 1
            k = SF.shape_key(gl)
            if k in seen_keys:
                stats["block_relabel_dup"] += 1
                continue
            seen_keys.add(k)
            stats["survive"] += 1
            shapes.append({"units": u, "potential": str(FR(pot4, 4)), "budget": 768 + (pot4 - 4 * (96 - 3 * N)) // 2,
                           "gadgets": [{"types": g["types"], "edges": g["edges"], "potential": g["potential"], "f": g["f"]} for g in gl]})
    print(f"slice {A}:{B} done: {json.dumps(dict(stats))}, {time.time() - t0:.0f} s", file=log, flush=True)
    shapes.sort(key=lambda s_: (-FR(s_["potential"]), len(s_["gadgets"]), sum(len(g["types"]) for g in s_["gadgets"])))
    json.dump({"N": N, "n": 72 + N, "bound": str(FR(bound4, 4)), "slice": [A, B], "stats": dict(stats),
               "wall_s": round(time.time() - t0, 1), "shapes": shapes}, open(out, "w"), indent=0)


def stage_merge(args):
    out = next(a.split("=")[1] for a in args if a.startswith("--out="))
    parts = [a for a in args if not a.startswith("--")]
    shapes, seen = [], set()
    stats = collections.Counter()
    meta = None
    for p in parts:
        d = json.load(open(p))
        meta = meta or {k: d[k] for k in ("N", "n", "bound")}
        for k, v in d["stats"].items():
            stats[k] += v
        for s in d["shapes"]:
            k = (s["units"], SF.shape_key(s["gadgets"]))
            if k in seen:
                stats["cross_slice_dup"] += 1
                continue
            seen.add(k)
            shapes.append(s)
    shapes.sort(key=lambda s: (-FR(s["potential"]), len(s["gadgets"]), sum(len(g["types"]) for g in s["gadgets"])))
    meta.update({"parts": parts, "stats": dict(stats), "shapes": shapes})
    json.dump(meta, open(out, "w"), indent=0)
    print(json.dumps({k: meta[k] for k in ("N", "n", "bound", "stats")}), len(shapes), "shapes", file=sys.stderr)


if __name__ == "__main__":
    {"enum": stage_enum, "match": stage_match, "merge": stage_merge}[sys.argv[1]](sys.argv[2:])
