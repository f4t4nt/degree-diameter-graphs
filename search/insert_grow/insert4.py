#!/usr/bin/env python3
"""Exact SAT completion: insert k new vertices into G' = G - Del with diameter <= D and max degree <= Delta.

Model (sound for any D; complete for paths with <= 2 new edges at one new vertex, one new-new edge chain, or
one old-old addition; exact for D = 2):
  y[j][s]  new vertex j adjacent to old vertex s (s in F = vertices of G' with a free slot)
  z[j][j'] new-new adjacency (--newadj free|indep|clique)
  w[s][s'] old-old addition between free vertices (--addable free, --maxadd M; default none)
  constraints: degree of every new vertex <= Delta (>= 1), free slots of every old vertex respected,
  every pair of the new graph at distance <= D (old-old far pairs of G' repaired, new-old, new-new).
Distances inside G' are precomputed (gtool dist / gtool far), so repair paths are enumerated explicitly.

usage: insert4.py --base G.edges --D 4 --delta 11 --k 1 [--delete "a-b,c-d"] [--delete-file F]
                  [--newadj free|indep|clique] [--addable none|free] [--maxadd M] [--time 600]
                  [--out sol.edges] [--maxslots S] [--mindeg m]
prints one summary line "RESULT ..." and, on a verified diameter-D solution, "!!! SAT n=<n> file=<out>".
"""
import argparse
import collections
import hashlib
import itertools
import os
import subprocess
import sys
import threading
import time

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Cadical153

HERE = os.path.dirname(os.path.abspath(__file__))
GTOOL = os.path.join(HERE, "gtool")
VERIFY = os.path.join(HERE, "verify_bfs")


def load(path):
    edges = set()
    n = 0
    for line in open(path):
        if line.strip() and not line.startswith("#"):
            a, b = map(int, line.split())
            edges.add((min(a, b), max(a, b)))
            n = max(n, a + 1, b + 1)
    return n, edges


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--D", type=int, default=4)
    ap.add_argument("--delta", type=int, required=True)
    ap.add_argument("--k", type=int, default=1)
    ap.add_argument("--delete", default="")
    ap.add_argument("--delete-file", default=None)
    ap.add_argument("--newadj", default="free")
    ap.add_argument("--addable", default="none")
    ap.add_argument("--maxadd", type=int, default=0)
    ap.add_argument("--time", type=float, default=600)
    ap.add_argument("--out", default=None)
    ap.add_argument("--mindeg", type=int, default=1)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    t0 = time.time()
    D, K, Delta = args.D, args.k, args.delta
    n, edges = load(args.base)
    dele = set()
    if args.delete:
        for tok in args.delete.split(","):
            a, b = map(int, tok.split("-"))
            dele.add((min(a, b), max(a, b)))
    if args.delete_file:
        for line in open(args.delete_file):
            if line.strip():
                a, b = map(int, line.split())
                dele.add((min(a, b), max(a, b)))
    missing = dele - edges
    if missing:
        sys.exit(f"deletions not in graph: {sorted(missing)[:5]}")
    gp = edges - dele
    adj = collections.defaultdict(set)
    for a, b in gp:
        adj[a].add(b)
        adj[b].add(a)
    slack = {v: Delta - len(adj[v]) for v in range(n) if Delta - len(adj[v]) > 0}
    F = sorted(slack)
    tag = args.tag or f"{os.path.basename(args.base)}_k{K}_del{len(dele)}"
    work = os.path.join(HERE, "runs", "sat_tmp")
    os.makedirs(work, exist_ok=True)
    gpath = os.path.join(work, f"{tag}.gp.edges")
    with open(gpath, "w") as f:
        for a, b in sorted(gp):
            f.write(f"{a} {b}\n")
    # far pairs of G'
    r = subprocess.run([GTOOL, "far", gpath, str(D)], capture_output=True, text=True)
    far = [tuple(map(int, l.split()[:2])) for l in r.stdout.splitlines() if l.strip()]
    # distances from F up to D-1
    fpath = os.path.join(work, f"{tag}.F.txt")
    with open(fpath, "w") as f:
        f.write("\n".join(map(str, F)) + "\n")
    r = subprocess.run([GTOOL, "dist", gpath, str(D - 1), fpath], capture_output=True, text=True)
    dF = {}          # dF[s] = dict u -> d (d <= D-1, u != s)
    nearF = collections.defaultdict(list)   # nearF[u] = [(s, d)] for s in F with d(s,u) <= D-1
    for line in r.stdout.splitlines():
        parts = list(map(int, line.split()))
        s = parts[0]
        dd = {}
        for i in range(1, len(parts), 2):
            dd[parts[i]] = parts[i + 1]
            nearF[parts[i]].append((s, parts[i + 1]))
        dF[s] = dd
    for s in F:
        nearF[s].append((s, 0))
    for tmp in (gpath, fpath):      # the G' copy and source list are only needed by gtool
        try:
            os.remove(tmp)
        except OSError:
            pass
    print(f"# {tag}: n={n} |G'|={len(gp)} deleted={len(dele)} far_pairs={len(far)} |F|={len(F)} slots={sum(slack.values())} ({time.time() - t0:.1f} s)", flush=True)
    if sum(slack.values()) < K * args.mindeg:
        print(f"RESULT {tag} UNSAT (not enough slots) 0.0s")
        return
    pool = IDPool()
    cnf = []
    y = {(j, s): pool.id(("y", j, s)) for j in range(K) for s in F}
    z = {}
    for j in range(K):
        for jj in range(j + 1, K):
            z[(j, jj)] = pool.id(("z", j, jj))
            if args.newadj == "indep":
                cnf.append([-z[(j, jj)]])
            elif args.newadj == "clique":
                cnf.append([z[(j, jj)]])
    def Z(a, b):
        return z[(min(a, b), max(a, b))]
    w = {}
    if args.addable == "free" and args.maxadd > 0:
        for i, s in enumerate(F):
            for s2 in F[i + 1:]:
                if s2 not in adj[s] and dF[s].get(s2, D) >= 2:
                    w[(s, s2)] = pool.id(("w", s, s2))
        cnf.extend(CardEnc.atmost(lits=list(w.values()), bound=args.maxadd, vpool=pool, encoding=EncType.seqcounter).clauses)
    def W(a, b):
        return w.get((min(a, b), max(a, b)))
    prod_cache = {}
    def prod(lits):
        """auxiliary p with p -> each literal (sound for positive use in disjunctions)"""
        key = tuple(sorted(lits))
        if len(key) == 1:
            return key[0]
        if key in prod_cache:
            return prod_cache[key]
        p = pool.id(("p", key))
        for l in key:
            cnf.append([-p, l])
        prod_cache[key] = p
        return p
    # degrees
    for j in range(K):
        lits = [y[(j, s)] for s in F] + [Z(j, jj) for jj in range(K) if jj != j]
        cnf.extend(CardEnc.atmost(lits=lits, bound=Delta, vpool=pool, encoding=EncType.seqcounter).clauses)
        if args.mindeg > 1:
            cnf.extend(CardEnc.atleast(lits=lits, bound=args.mindeg, vpool=pool, encoding=EncType.seqcounter).clauses)
        else:
            cnf.append(lits)
    for s in F:
        lits = [y[(j, s)] for j in range(K)] + [v for (a, b), v in w.items() if s in (a, b)]
        if len(lits) > slack[s]:
            cnf.extend(CardEnc.atmost(lits=lits, bound=slack[s], vpool=pool, encoding=EncType.seqcounter).clauses)
    # symmetry breaking for interchangeable new vertices (newadj indep/clique: all new vertices are alike):
    # new vertex j's smallest F-neighbour index is < that of j+1  (encoded via prefix variables)
    if K > 1 and args.newadj in ("indep", "clique"):
        pre = {}
        for j in range(K):
            acc = None
            for i, s in enumerate(F):
                v = pool.id(("pre", j, i))   # v <-> OR_{i' <= i} y[j][F[i']]
                if acc is None:
                    cnf.append([-v, y[(j, s)]]); cnf.append([v, -y[(j, s)]])
                else:
                    cnf.append([-v, acc, y[(j, s)]]); cnf.append([v, -acc]); cnf.append([v, -y[(j, s)]])
                pre[(j, i)] = v
                acc = v
        for j in range(K - 1):
            for i in range(len(F)):
                # if j+1 has a neighbour among F[0..i] then j has one among F[0..i-1]
                if i == 0:
                    cnf.append([-y[(j + 1, F[0])]])
                else:
                    cnf.append([-pre[(j + 1, i)], pre[(j, i - 1)]])
    # (c) old-old far pairs of G'
    nrep = 0
    for (u, v) in far:
        cl = []
        nu = nearF[u]   # (s, d(s,u))
        nv = nearF[v]
        for j in range(K):
            for (s, a) in nu:
                if a > D - 2:
                    continue
                for (s2, b) in nv:
                    if s2 != s and a + b <= D - 2:
                        cl.append(prod([y[(j, s)], y[(j, s2)]]))
            for jj in range(K):
                if jj == j:
                    continue
                for (s, a) in nu:
                    if a > D - 3:
                        continue
                    for (s2, b) in nv:
                        if a + b <= D - 3:
                            cl.append(prod([y[(j, s)], Z(j, jj), y[(jj, s2)]]))
            if D >= 4 and u in slack and v in slack:
                for jj in range(K):
                    for jjj in range(K):
                        if len({j, jj, jjj}) == 3:
                            cl.append(prod([y[(j, u)], Z(j, jj), Z(jj, jjj), y[(jjj, v)]]))
        if w:
            for (s, a) in nu:
                for (s2, b) in nv:
                    if s != s2 and a + 1 + b <= D and W(s, s2) is not None:
                        cl.append(W(s, s2))
            if D == 2 and u in slack and v in slack:
                for s in F:
                    if W(u, s) is not None and W(s, v) is not None:
                        cl.append(prod([W(u, s), W(s, v)]))
            # mixed: u -a- s - x_j - s2 - s3 -b- v with a + b <= D - 3
            for j in range(K):
                for (s, a) in nu:
                    if a > D - 3:
                        continue
                    for (s3, b) in nv:
                        if a + b > D - 3:
                            continue
                        for s2 in F:
                            if s2 != s and W(s2, s3) is not None:
                                cl.append(prod([y[(j, s)], y[(j, s2)], W(s2, s3)]))
        if not cl:
            print(f"RESULT {tag} UNSAT (far pair {u}-{v} unrepairable) {time.time() - t0:.1f}s")
            return
        cnf.append(cl)
        nrep += len(cl)
    # (d) new-old pairs
    for j in range(K):
        for u in range(n):
            cl = [y[(j, s)] for (s, a) in nearF[u] if a <= D - 1]
            for jj in range(K):
                if jj == j:
                    continue
                for (s, a) in nearF[u]:
                    if a <= D - 2:
                        cl.append(prod([Z(j, jj), y[(jj, s)]]))
                for jjj in range(K):
                    if jjj in (j, jj):
                        continue
                    for (s, a) in nearF[u]:
                        if a <= D - 3:
                            cl.append(prod([Z(j, jj), Z(jj, jjj), y[(jjj, s)]]))
            if w:
                for (s2, a) in nearF[u]:
                    if a <= D - 2:
                        for s in F:
                            if s != s2 and W(s, s2) is not None:
                                cl.append(prod([y[(j, s)], W(s, s2)]))
            if not cl:
                print(f"RESULT {tag} UNSAT (vertex {u} unreachable from new vertex) {time.time() - t0:.1f}s")
                return
            cnf.append(cl)
    # (e) new-new pairs
    for j in range(K):
        for jj in range(j + 1, K):
            cl = [Z(j, jj)]
            for jjj in range(K):
                if jjj not in (j, jj):
                    cl.append(prod([Z(j, jjj), Z(jjj, jj)]))
            for s in F:
                for (s2, b) in nearF[s]:
                    if b <= D - 2 and s2 in slack:
                        cl.append(prod([y[(j, s)], y[(jj, s2)]]))
            cnf.append(cl)
    nv_ = pool.top
    print(f"# {tag}: vars={nv_} clauses={len(cnf)} repair_lits={nrep} build {time.time() - t0:.1f} s", flush=True)
    solver = Cadical153(bootstrap_with=cnf)
    # time cap by repeated conflict budgets (solver.interrupt() is not honoured reliably by this CaDiCaL build)
    t1 = time.time()
    res = None
    budget = 20000
    while True:
        solver.conf_budget(budget)
        res = solver.solve_limited()
        if res is not None or time.time() - t1 > args.time:
            break
        budget = min(budget * 2, 500000)
    status = {True: "SAT", False: "UNSAT", None: "TIMEOUT"}[res]
    print(f"RESULT {tag} {status} solve {time.time() - t1:.1f}s total {time.time() - t0:.1f}s vars={nv_} clauses={len(cnf)} far={len(far)} |F|={len(F)}", flush=True)
    if not res:
        return
    model = set(l for l in solver.get_model() if l > 0)
    new_edges = set(gp)
    for (j, s), v in y.items():
        if v in model:
            new_edges.add((s, n + j))
    for (j, jj), v in z.items():
        if v in model:
            new_edges.add((n + j, n + jj))
    for (s, s2), v in w.items():
        if v in model:
            new_edges.add((s, s2))
    out = args.out or os.path.join(HERE, "runs", f"sol_{tag}_{n + K}.edges")
    with open(out, "w") as f:
        for a, b in sorted(new_edges):
            f.write(f"{a} {b}\n")
    ver = subprocess.run([VERIFY, out, str(D)], capture_output=True, text=True).stdout
    print(ver)
    sha = hashlib.sha256(open(out, "rb").read()).hexdigest()
    ok = f"unordered pairs at distance > {D}: 0" in ver and f"max degree: " in ver and int(ver.split("max degree: ")[1].split()[0]) <= Delta
    for j in range(K):
        print(f"# new vertex {n + j}: N = {sorted(b if a == n + j else a for a, b in new_edges if n + j in (a, b))}")
    if ok:
        print(f"!!! SAT n={n + K} file={out} sha256={sha}", flush=True)
    else:
        print(f"!!! MODEL FAILED VERIFICATION file={out}", flush=True)


if __name__ == "__main__":
    main()
