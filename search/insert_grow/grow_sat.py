#!/usr/bin/env python3
"""SAT-assisted iterated insertion for the (11,4) W(9)-double family (continues grow.c's output).

Each step deletes the two point-twin (or point-spread-line) edges at the two copies of a base point p and one
line-line or twin-twin edge of a base line L (candidates ordered by total crit_4, allowing crit_4 <= CMAX), and
asks insert4.py for ONE new vertex using those six endpoints and every free slot of the graph (earlier new
vertices); the far pairs created (<= 3 * CMAX) must be repaired exactly.  First SAT is accepted; the step is
re-verified by verify_bfs through insert4.  usage: grow_sat.py START.edges OUTPREFIX [max_steps] [cmax] [max_trials]"""
import collections
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
NP, NL, NT = 1640, 1640, 164


def typ(v):
    return "P" if v < NP else ("L" if v < NP + NL else ("T" if v < NP + NL + NT else "X"))


def crit_of(path):
    out = subprocess.run([os.path.join(HERE, "gtool"), "crit", path, "4"], capture_output=True, text=True).stdout
    crit = {}
    for line in out.splitlines():
        if line.startswith("#"):
            continue
        u, v, c1, c2, c3, c4 = map(int, line.split())
        crit[(u, v)] = c4
    return crit


def main():
    start, prefix = sys.argv[1], sys.argv[2]
    max_steps = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
    cmax = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    max_trials = int(sys.argv[5]) if len(sys.argv) > 5 else 60
    cur = f"{prefix}_cur.edges"
    subprocess.run(["cp", start, cur])
    log = open(f"{prefix}.log", "a")
    t0 = time.time()
    for step in range(max_steps):
        edges = set()
        adj = collections.defaultdict(set)
        n = 0
        for line in open(cur):
            a, b = map(int, line.split())
            edges.add((a, b))
            adj[a].add(b)
            adj[b].add(a)
            n = max(n, a + 1, b + 1)
        crit = crit_of(cur)
        # point candidates: for each base point i, the cheapest deletable edge at each copy (to a twin or to a line
        # whose other endpoint is a line copy with crit <= cmax)
        pts = []
        for i in range(820):
            tot, chosen = 0, []
            for v in (2 * i, 2 * i + 1):
                best = None
                for w in adj[v]:
                    if typ(w) not in ("T", "L"):
                        continue
                    e = (min(v, w), max(v, w))
                    c = crit[e]
                    if c <= cmax and (best is None or c < best[0]):
                        best = (c, w)
                if best is None:
                    break
                tot += best[0]
                chosen.append(best[1])
            if len(chosen) == 2:
                pts.append((tot, i, chosen))
        pts.sort()
        lns = []
        for j in range(820):
            L0, L1 = NP + 2 * j, NP + 2 * j + 1
            if (L0, L1) in edges and crit[(L0, L1)] <= cmax:
                lns.append((crit[(L0, L1)], L0, L1))
        for j in range(82):
            T0, T1 = NP + NL + 2 * j, NP + NL + 2 * j + 1
            if (T0, T1) in edges and crit[(T0, T1)] <= cmax:
                lns.append((crit[(T0, T1)], T0, T1))
        lns.sort()
        print(f"step {step + 1}: n={n} point candidates {len(pts)} (best total crit {pts[0][0] if pts else None}), line/twin candidates {len(lns)} (best {lns[0][0] if lns else None}); {time.time() - t0:.0f} s", flush=True)
        if not pts or not lns:
            break
        trials = 0
        accepted = False
        # interleave: cheapest combinations first (sum of crit)
        combos = sorted(((pc + lc, pi, pch, l0, l1) for (pc, pi, pch) in pts[:40] for (lc, l0, l1) in lns[:40]), key=lambda t: t[0])
        for (tot, i, chosen, L0, L1) in combos:
            if trials >= max_trials:
                break
            a, b = 2 * i, 2 * i + 1
            if L0 in adj[a] or L0 in adj[b] or L1 in adj[a] or L1 in adj[b]:
                continue
            dele = f"{a}-{chosen[0]},{b}-{chosen[1]},{L0}-{L1}"
            tag = f"gs_{n}_{i}_{L0}"
            out = f"{prefix}_trial.edges"
            r = subprocess.run([PY, os.path.join(HERE, "insert4.py"), "--base", cur, "--D", "4", "--delta", "11", "--k", "1",
                                "--delete", dele, "--time", "60", "--tag", tag, "--out", out], capture_output=True, text=True)
            trials += 1
            res = [l for l in r.stdout.splitlines() if l.startswith("RESULT")]
            if any(l.startswith("!!! SAT") for l in r.stdout.splitlines()):
                nb = [l for l in r.stdout.splitlines() if l.startswith("# new vertex")]
                msg = f"step {step + 1}: n={n + 1} point {i} line/twin {L0}-{L1} deleted {dele} (crit total {tot}); {nb[0] if nb else ''}; trial {trials}"
                print(msg, flush=True)
                log.write(msg + "\n")
                log.flush()
                subprocess.run(["cp", out, cur])
                subprocess.run(["cp", out, f"{prefix}_{n + 1}.edges"])
                if os.path.exists(f"{prefix}_{n}.edges"):
                    os.remove(f"{prefix}_{n}.edges")
                accepted = True
                break
            elif trials <= 2:
                print("   ", res[0][:160] if res else r.stdout[-300:], flush=True)
        if not accepted:
            print(f"step {step + 1}: no SAT among {trials} trials", flush=True)
            break
    print("final", cur)


if __name__ == "__main__":
    main()
