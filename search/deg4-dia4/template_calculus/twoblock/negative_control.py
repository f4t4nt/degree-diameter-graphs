"""Negative control for twoblock_sat.py: perturb a known 104 inside W (a degree-preserving double
edge swap among plain edges, or a swap of two port attachments), check by BFS that the diameter
grows, and verify that the encoder rejects the perturbed graph under full assumptions (UNSAT)."""
import collections
import json
import random
import sys

import twoblock_sat as TS


def diameter(edges, n):
    return TS.bfs_check(edges, n)


def main():
    path = sys.argv[1]
    edges = TS.relabel_known(path)
    n = len({v for e in edges for v in e})
    es = set((min(u, v), max(u, v)) for u, v in edges)
    rng = random.Random(1)
    W = list(range(48, n))
    out = []
    trials = 0
    while len(out) < 3 and trials < 200:
        trials += 1
        plain = [e for e in es if e[0] >= 48 and e[1] >= 48]
        (a, b), (c, d) = rng.sample(plain, 2)
        if len({a, b, c, d}) < 4:
            continue
        new = [(min(a, c), max(a, c)), (min(b, d), max(b, d))]
        if any(e in es for e in new):
            continue
        es2 = (es - {(a, b), (c, d)}) | set(new)
        chk = diameter(sorted(es2), n)
        if chk.get("pairs_beyond_4", 0) == 0:
            continue
        model = TS.Model(n, "general")
        res = model.solve(assumptions=model.assumptions_from_graph(sorted(es2)))
        out.append({"swap": [[a, b], [c, d]], "bfs": chk, "sat_under_assumptions": res})
        print(out[-1], file=sys.stderr, flush=True)
    json.dump(out, open("runs/negative_control.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
