#!/usr/bin/env python3
"""Self-contained rebuild of a grown graph from the base edge list and the step list (no search code).

usage: rebuild_chain.py BASE.edges OUT.edges STEPFILE...
STEPFILE: grow.c logs ("step k: new vertex x (degree d) deleted a-b c-d ; N(x) = s1 s2 ... ; trial t") and grow_sel /
q23_job witness files ("# step k: STEP n=<n> new vertex x (degree d) deleted a-b ... ; N(x) = s1 ... ; attempt a ; t s"),
in chain order.  Each step: delete the listed edges (they must exist), add vertex x = current n with the listed
neighbours (x must equal n).  Writes OUT.edges as sorted "u v" lines (u < v) and prints n, m and the sha256.
"""
import hashlib
import re
import sys


def main():
    base, out = sys.argv[1], sys.argv[2]
    edges = set()
    n = 0
    for line in open(base):
        if line.strip():
            a, b = map(int, line.split())
            edges.add((min(a, b), max(a, b)))
            n = max(n, a + 1, b + 1)
    steps = 0
    for path in sys.argv[3:]:
        for line in open(path):
            m = re.search(r"new vertex (\d+) \(degree (\d+)\) deleted (.*?) ; N\(x\) = (.*?) ;", line)
            if not m:
                continue
            x, dg = int(m.group(1)), int(m.group(2))
            dels = [tuple(map(int, t.split("-"))) for t in m.group(3).split()]
            nbrs = list(map(int, m.group(4).split()))
            if x != n:
                sys.exit(f"{path}: step for vertex {x} but current order is {n}")
            for a, b in dels:
                e = (min(a, b), max(a, b))
                if e not in edges:
                    sys.exit(f"{path}: deleted edge {e} not present at step of vertex {x}")
                edges.remove(e)
            if len(nbrs) != dg or len(set(nbrs)) != dg:
                sys.exit(f"{path}: neighbourhood of {x} has {len(nbrs)} entries, degree {dg}")
            for s in nbrs:
                if s >= n:
                    sys.exit(f"{path}: neighbour {s} of {x} does not exist yet")
                edges.add((s, x))
            n += 1
            steps += 1
    with open(out, "w") as f:
        for a, b in sorted(edges):
            f.write(f"{a} {b}\n")
    sha = hashlib.sha256(open(out, "rb").read()).hexdigest()
    print(f"rebuilt {steps} steps: n={n} m={len(edges)} sha256={sha} -> {out}")


if __name__ == "__main__":
    main()
