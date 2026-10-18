#!/usr/bin/env python3
"""Iterated insertion with the DELETION SET CHOSEN BY SAT ("selection model"), one new vertex per step.

Step on the current graph G (max degree Delta, diameter D):
  pool P = the edges with crit_D <= cmax (at most `pool` of them, lowest crit first); for every e in P the exact set
  F_e of unordered pairs that deleting e alone pushes beyond D (farsets.c: candidates = pairs all of whose geodesics
  use e, re-checked by BFS in G - e).  S = endpoints of P + vertices with a free slot.
  Variables: del_e (e in P), y_s (s in S: the new vertex x is adjacent to s), reach r_{w,t} (w within t of N(x)).
  Clauses: y_s needs a slot (free slot or a deleted pool edge at s); sum y_s <= Delta; for every e and (u, v) in
  F_e: del_e -> (r_{u,0} & r_{v,D-2}) | (r_{u,1} & r_{v,D-3}) | ... (a repair path u .. s - x - s' .. v of length
  <= D); r_{w,t} -> OR of y_s with d_G(s, w) <= t; cover: every vertex w has r_{w,D-1}.
  Distances are G-distances (a route might use a deleted edge), and far pairs created jointly by several deletions
  are not modelled, so every model is verified with verify_bfs; a failed model is blocked and the solver re-run.
usage: grow_sel.py START.edges OUTPREFIX --delta 11 --D 4 [--cmax 4] [--pool 200] [--fcap 60] [--steps 1000]
                   [--time 120] [--blocks 20] [--budget 3300] [--threshold 3511] [--seed 1] [--maxdel 0]
prints "STEP n=<n> ..." per accepted step and one final "FOUND n=<n> deg=<Delta> D=<D> file=<edges> sha256=<sha>".
"""
import argparse
import collections
import hashlib
import os
import random
import shutil
import subprocess
import sys
import time

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Cadical153

HERE = os.path.dirname(os.path.abspath(__file__))
GTOOL = os.path.join(HERE, "gtool")
FARSETS = os.path.join(HERE, "farsets")
VERIFY = os.path.join(HERE, "verify_bfs")


def load(path):
    adj = collections.defaultdict(set)
    n = 0
    for line in open(path):
        if line.strip():
            a, b = map(int, line.split())
            adj[a].add(b)
            adj[b].add(a)
            n = max(n, a + 1, b + 1)
    return n, adj


def crit_of(path, D):
    out = subprocess.run([GTOOL, "crit", path, str(D)], capture_output=True, text=True).stdout
    crit = {}
    for line in out.splitlines():
        if line.startswith("#"):
            continue
        p = list(map(int, line.split()))
        crit[(p[0], p[1])] = p[-1]
    return crit


def dists(path, R, sources):
    """{s: {w: d}} for d <= R (w != s), via gtool dist"""
    spath = path + ".src"
    with open(spath, "w") as f:
        f.write("\n".join(map(str, sources)) + "\n")
    out = subprocess.run([GTOOL, "dist", path, str(R), spath], capture_output=True, text=True).stdout
    os.remove(spath)
    res = {}
    for line in out.splitlines():
        p = list(map(int, line.split()))
        res[p[0]] = {p[i]: p[i + 1] for i in range(1, len(p), 2)}
    return res


def farsets(path, D, pool):
    ppath = path + ".pool"
    with open(ppath, "w") as f:
        f.write("".join(f"{a} {b}\n" for a, b in pool))
    out = subprocess.run([FARSETS, path, str(D), ppath], capture_output=True, text=True).stdout
    os.remove(ppath)
    res = {}
    for line in out.splitlines():
        head, _, tail = line.partition(":")
        a, b = map(int, head.split())
        p = list(map(int, tail.split()))
        res[(a, b)] = [(p[i], p[i + 1]) for i in range(0, len(p), 2)]
    return res


def verify(path, D):
    out = subprocess.run([VERIFY, path, str(D)], capture_output=True, text=True).stdout
    far = int(out.split(f"distance > {D}:")[1].split()[0])
    maxd = int(out.split("max degree:")[1].split()[0])
    order = int(out.split("order:")[1].split()[0])
    dup = int(out.split("duplicate edges:")[1].split()[0])
    return out, order, maxd, far, dup


def step(cur, D, Delta, args, rng, log, cmax, poolsz, fcap, maxxx):
    t0 = time.time()
    n, adj = load(cur)
    crit = crit_of(cur, D)
    cand = sorted(((c, u, v) for (u, v), c in crit.items() if c <= cmax))
    if args.seed > 1:
        # diversification between seeds: shuffle within equal crit (truncated pools differ), drop a random fraction of
        # the pool, and (below) shuffle the variable creation and clause order seen by the solver
        rng.shuffle(cand)
        cand.sort(key=lambda t: t[0])
        cand = [t for t in cand if rng.random() >= args.drop]
    pool = [(u, v) for c, u, v in cand[:poolsz]]
    if args.seed > 1:
        rng.shuffle(pool)
    F = farsets(cur, D, pool)
    pool = [e for e in pool if len(F[e]) <= fcap]
    slack = {v: Delta - len(adj[v]) for v in range(n) if len(adj[v]) < Delta}
    S = sorted({u for u, v in pool} | {v for u, v in pool} | set(slack))
    if args.seed > 1:
        rng.shuffle(S)
    dS = dists(cur, D - 1, S)
    # near[w] = [(s, d)] for s in S with d(s, w) <= D-1
    near = collections.defaultdict(list)
    for s in S:
        near[s].append((s, 0))
        for w, d in dS[s].items():
            near[w].append((s, d))
    nfar = sum(len(F[e]) for e in pool)
    log(f"  n={n} pool {len(pool)} edges (crit<={cmax}, |F_e|<={fcap}, {nfar} far pairs), |S|={len(S)}, free slots {sum(slack.values())} ({time.time() - t0:.1f} s)")
    pool_at = collections.defaultdict(list)
    for e in pool:
        pool_at[e[0]].append(e)
        pool_at[e[1]].append(e)
    vp = IDPool()
    cnf = []
    dl = {e: vp.id(("del", e)) for e in pool}
    y = {s: vp.id(("y", s)) for s in S}
    reach = {}

    def r(w, t):
        """variable implying: some neighbour s of x has d(s, w) <= t"""
        key = (w, t)
        if key in reach:
            return reach[key]
        lits = [y[s] for s, d in near[w] if d <= t]
        if not lits:
            reach[key] = None
            return None
        v = vp.id(("r", w, t))
        cnf.append([-v] + lits)
        reach[key] = v
        return v

    prods = {}

    def prod(a, b):
        key = (min(a, b), max(a, b))
        if key in prods:
            return prods[key]
        p = vp.id(("p", key))
        cnf.append([-p, a])
        cnf.append([-p, b])
        prods[key] = p
        return p

    # slots; a deletion must free a slot that x uses (no gratuitous deletions: they only add joint far pairs)
    for s in S:
        if slack.get(s, 0) == 0:
            cnf.append([-y[s]] + [dl[e] for e in pool_at[s]])
    for e in pool:
        cnf.append([-dl[e], y[e[0]], y[e[1]]])
    cnf.extend(CardEnc.atmost(lits=list(y.values()), bound=Delta, vpool=vp, encoding=EncType.seqcounter).clauses)
    cnf.append(list(y.values()))
    if args.maxdel:
        cnf.extend(CardEnc.atmost(lits=list(dl.values()), bound=args.maxdel, vpool=vp, encoding=EncType.seqcounter).clauses)
    if maxxx >= 0 and args.n0:
        # cap on deletions of edges between NEW vertices (index >= n0): they are the layer's own connectivity
        xx = [dl[e] for e in pool if e[0] >= args.n0 and e[1] >= args.n0]
        if maxxx == 0:
            cnf.extend([-v] for v in xx)
        elif len(xx) > maxxx:
            cnf.extend(CardEnc.atmost(lits=xx, bound=maxxx, vpool=vp, encoding=EncType.seqcounter).clauses)
    # repairs: reach variables per (deleted edge e, endpoint u, radius t) from BFS in G - e (exact for a single
    # deletion; a repair route u .. s - x - s' .. v of length <= D with s = s' would make d_{G-e}(u, v) <= D-2)
    Sset = set(S)
    forbidden = 0
    for e in pool:
        ends = {u for u, v in F[e]} | {v for u, v in F[e]}
        rad = {}   # (u, t) -> list of y literals
        for u in ends:
            dist = {u: 0}
            frontier = [u]
            for t in range(1, D - 1):
                nxt = []
                for x in frontier:
                    for w in adj[x]:
                        if (x, w) in (e, e[::-1]) or w in dist:
                            continue
                        dist[w] = t
                        nxt.append(w)
                frontier = nxt
            for t in range(0, D - 1):
                rad[(u, t)] = [y[s] for s, d in dist.items() if s in Sset and d <= t]
        rv_ = {}

        def re_(u, t):
            if (u, t) in rv_:
                return rv_[(u, t)]
            lits = rad[(u, t)]
            if not lits:
                rv_[(u, t)] = None
                return None
            if len(lits) == 1:
                rv_[(u, t)] = lits[0]
                return lits[0]
            v_ = vp.id(("re", e, u, t))
            cnf.append([-v_] + lits)
            rv_[(u, t)] = v_
            return v_

        for (u, v) in F[e]:
            cl = [-dl[e]]
            for t in range(0, D - 1):
                ru, rv = re_(u, t), re_(v, D - 2 - t)
                if ru is not None and rv is not None:
                    cl.append(prod(ru, rv))
            if len(cl) == 1:
                forbidden += 1
            cnf.append(cl)
    # cover
    for w in range(n):
        rw = r(w, D - 1)
        if rw is None:
            log(f"  vertex {w} is not within {D - 1} of any candidate neighbour: UNSAT")
            return None
        cnf.append([rw])
    log(f"  model: {vp.top} vars, {len(cnf)} clauses, {forbidden} pool edges undeletable ({time.time() - t0:.1f} s)")
    if args.seed > 1:
        rng.shuffle(cnf)
    solver = Cadical153(bootstrap_with=cnf)
    t1 = time.time()
    for attempt in range(args.blocks + 1):
        budget = 20000
        res = None
        while True:
            solver.conf_budget(budget)
            res = solver.solve_limited()
            if res is not None or time.time() - t1 > args.time:
                break
            budget = min(budget * 2, 500000)
        if res is None:
            log(f"  TIMEOUT after {time.time() - t1:.0f} s (attempt {attempt})")
            return None
        if res is False:
            log(f"  UNSAT ({time.time() - t1:.1f} s, attempt {attempt}): no single-vertex insertion over this pool")
            return None
        model = set(l for l in solver.get_model() if l > 0)
        dels = [e for e in pool if dl[e] in model]
        nbrs = [s for s in S if y[s] in model]
        edges = set()
        for a in adj:
            for b in adj[a]:
                if a < b:
                    edges.add((a, b))
        for e in dels:
            edges.discard(e)
        for s in nbrs:
            edges.add((s, n))
        out = f"{args.prefix}_trial.edges"
        with open(out, "w") as f:
            for a, b in sorted(edges):
                f.write(f"{a} {b}\n")
        ver, order_, maxd, far, dup = verify(out, D)
        if far == 0 and dup == 0 and maxd <= Delta and order_ == n + 1:
            sha = hashlib.sha256(open(out, "rb").read()).hexdigest()
            return out, dels, nbrs, ver, sha, attempt
        log(f"  model failed verification (far {far}, maxdeg {maxd}, dup {dup}; {len(dels)} deletions, degree {len(nbrs)}); blocking the deletion set (attempt {attempt})")
        solver.add_clause([-dl[e] for e in dels])
    log(f"  gave up after {args.blocks} blocked models")
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("start")
    ap.add_argument("prefix")
    ap.add_argument("--delta", type=int, required=True)
    ap.add_argument("--D", type=int, required=True)
    ap.add_argument("--cmax", type=int, default=4)
    ap.add_argument("--pool", type=int, default=200)
    ap.add_argument("--fcap", type=int, default=60)
    ap.add_argument("--steps", type=int, default=1000)
    ap.add_argument("--time", type=float, default=120)
    ap.add_argument("--blocks", type=int, default=20)
    ap.add_argument("--budget", type=float, default=3300)
    ap.add_argument("--threshold", type=int, default=0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--maxdel", type=int, default=0)
    ap.add_argument("--keep-last", action="store_true", help="keep only the latest PREFIX_<n>.edges (the witness file records every step)")
    ap.add_argument("--drop", type=float, default=0.15, help="fraction of the pool dropped at random per step when seed > 1")
    ap.add_argument("--levels", type=int, default=5, help="escalation levels tried after an UNSAT / timeout (1 = none)")
    ap.add_argument("--restarts", type=int, default=0, help="after the whole ladder fails, retry with a fresh random pool drop (seed > 1 only), this many times")
    ap.add_argument("--n0", type=int, default=0, help="vertices >= n0 are the inserted layer; with n0 > 0 the ladder first caps new-new deletions at 0, then 2")
    args = ap.parse_args()
    rng = random.Random(args.seed)
    Delta, D = args.delta, args.D
    t0 = time.time()
    cur = f"{args.prefix}_cur.edges"
    shutil.copy(args.start, cur)
    wit = open(f"{args.prefix}.witness.txt", "a")
    base_sha = hashlib.sha256(open(args.start, "rb").read()).hexdigest()
    wit.write(f"# grow_sel start {args.start} sha256 {base_sha} delta {Delta} D {D} cmax {args.cmax} pool {args.pool} fcap {args.fcap} seed {args.seed}\n")
    wit.flush()

    def log(msg):
        print(msg, flush=True)

    best = None
    # escalation ladder after an UNSAT / timeout: (cmax and fcap multiplier, pool multiplier, cap on new-new deletions;
    # -1 = no cap); with --n0 0 the caps are inactive and the ladder is the plain pool escalation
    LADDER = [(1, 1, 0), (1, 1, 2), (1, 1, -1), (2, 1.5, -1), (4, 2, -1)] if args.n0 else [(1, 1, -1), (2, 1.5, -1), (4, 2, -1)]
    if args.levels < len(LADDER):
        LADDER = LADDER[-args.levels:] if args.n0 == 0 else LADDER[:args.levels]
    level = 0
    restarts_left = args.restarts
    for st in range(args.steps):
        if time.time() - t0 > args.budget:
            log(f"budget exhausted before step {st + 1}")
            break
        mult, pmult, maxxx = LADDER[level]
        log(f"step {st + 1} ({time.time() - t0:.0f} s) level {level} (cmax x{mult}, maxxx {maxxx})")
        res = step(cur, D, Delta, args, rng, log, args.cmax * mult, int(args.pool * pmult), args.fcap * mult, maxxx)
        if res is None:
            if level >= len(LADDER) - 1:
                if restarts_left <= 0 or args.seed <= 1:
                    break
                restarts_left -= 1
                rng.seed(args.seed * 1000 + restarts_left)   # a fresh random pool drop / shuffle, same ladder
                log(f"  ladder exhausted: restart with a fresh random pool ({restarts_left} restarts left)")
                level = 0
                continue
            level += 1
            continue
        level = 0
        out, dels, nbrs, ver, sha, attempt = res
        n1 = int(ver.split("order:")[1].split()[0])
        dst = f"{args.prefix}_{n1}.edges"
        shutil.copy(out, dst)
        shutil.copy(out, cur)
        dele = " ".join(f"{a}-{b}" for a, b in dels)
        msg = f"STEP n={n1} new vertex {n1 - 1} (degree {len(nbrs)}) deleted {dele} ; N(x) = {' '.join(map(str, nbrs))} ; attempt {attempt} ; {time.time() - t0:.0f} s"
        log(msg)
        wit.write(f"# step {st + 1}: {msg}\n# file {dst} sha256 {sha}\n")
        wit.write("".join(f"#   {l}\n" for l in ver.strip().splitlines()))
        wit.flush()
        if args.keep_last and best:
            try:
                os.remove(best[1])
            except OSError:
                pass
        best = (n1, dst, sha)
    if best and best[0] > args.threshold:
        print(f"FOUND n={best[0]} deg={Delta} D={D} file={os.path.abspath(best[1])} sha256={best[2]}", flush=True)
    print(f"final n={best[0] if best else 'start'} ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
