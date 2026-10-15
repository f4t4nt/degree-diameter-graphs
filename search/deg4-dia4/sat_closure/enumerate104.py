"""Enumerate every diameter-4 realisation of a skeleton at n = 104 (all split assignments), reduce to
isomorphism classes, and check which of the classes in ../classes/ are among them (end-to-end
validation of the tightness theory, the skeleton enumeration and the SAT encoder).
Usage: python enumerate104.py <skeleton index> [phase1|phase2|both]

Phase 1 streams labelled solutions (option index per unit) to runs/all104_skel<k>_solutions.txt.
Phase 2 reduces them by Weisfeiler-Lehman hash plus exact isomorphism on hash collisions, keeping
one representative graph per class in memory."""
import collections
import json
import pathlib
import sys
import time

import networkx as nx
from pysat.solvers import Cadical153

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "template_calculus"))
import sat105  # noqa: E402
import verify  # noqa: E402


def phase1(skel, enc, path):
    solver = Cadical153(bootstrap_with=enc.cnf.clauses)
    sel_vars = [(u, o, v) for u in enc.sel for o, v in enumerate(enc.sel[u])]
    count = 0
    t0 = time.time()
    with open(path, "w") as f:
        while solver.solve():
            model = solver.get_model()
            pos = set(l for l in model if l > 0)
            chosen = [(u, o, v) for (u, o, v) in sel_vars if v in pos]
            solver.add_clause([-v for (_, _, v) in chosen])
            f.write(" ".join(str(o) for (_, o, _) in sorted(chosen)) + "\n")
            count += 1
            if count % 2000 == 0:
                print(count, "solutions", round(time.time() - t0, 1), "s", file=sys.stderr, flush=True)
    solver.delete()
    print("labelled solutions:", count, round(time.time() - t0, 1), "s", file=sys.stderr, flush=True)
    return count


def edges_from_options(enc, opts):
    edges = set()
    for v, s in enc.fixed.items():
        for w in s:
            edges.add(tuple(sorted((str(v), str(w)))))
    for u, o in enumerate(opts):
        kind, holds, extra, active = enc.options[u][o]
        for k in range(3):
            for t in holds[k]:
                edges.add(tuple(sorted((str(("C", u, k)), str(t)))))
        edges.add(tuple(sorted((str(("C", u, 0)), str(("C", u, 1))))))
        for (i, j) in extra:
            edges.add(tuple(sorted((str(("C", u, i)), str(("C", u, j))))))
    return sorted(edges)


def phase2(skel, enc, path):
    buckets = {}
    n = 0
    t0 = time.time()
    for line in open(path):
        opts = [int(x) for x in line.split()]
        edges = edges_from_options(enc, opts)
        G = nx.Graph()
        G.add_edges_from(edges)
        h = nx.weisfeiler_lehman_graph_hash(G, iterations=4)
        reps = buckets.setdefault(h, [])
        for rep in reps:
            if nx.is_isomorphic(G, rep[0]):
                rep[1] += 1
                break
        else:
            v = sat105.verify_edges(edges)
            assert v["diameter"] == 4 and v["bad_pairs"] == 0, v
            reps.append([G, 1])
        n += 1
        if n % 5000 == 0:
            print(n, "processed,", sum(len(r) for r in buckets.values()), "classes", round(time.time() - t0, 1), "s",
                  file=sys.stderr, flush=True)
    classes = [(G, cnt) for reps in buckets.values() for G, cnt in reps]
    print("isomorphism classes among solutions:", len(classes), file=sys.stderr, flush=True)
    corpus = verify.load_corpus()
    nk = len(corpus["representatives"])
    missing = []
    for i, rep in enumerate(corpus["representatives"]):
        Gk = nx.Graph()
        Gk.add_edges_from([tuple(e) for e in rep["edges_sorted"]])
        hk = nx.weisfeiler_lehman_graph_hash(Gk, iterations=4)
        if not any(nx.is_isomorphic(G, Gk) for G, _ in buckets.get(hk, [])):
            missing.append(i)
    print("known classes recovered:", nk - len(missing), f"of {nk}; missing:", missing, file=sys.stderr, flush=True)
    reps_edges = [sorted(tuple(sorted(e)) for e in G.edges()) for G, _ in classes]
    out = {"skeleton": skel, "labelled_solutions": n, "isomorphism_classes": len(classes),
           "known_recovered": nk - len(missing), "known_missing": missing,
           "class_sizes": sorted([cnt for _, cnt in classes], reverse=True)}
    with open(f"runs/all104_skel{skel}_classes.json", "w") as f:
        json.dump([{"edges": e} for e in reps_edges], f)
    return out


def main():
    skel = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    units = sat105.load_skeletons()[skel]
    enc = sat105.Encoder(units, "104")
    enc.encode_distances()
    path = f"runs/all104_skel{skel}_solutions.txt"
    mode = sys.argv[2] if len(sys.argv) > 2 else "both"
    if mode != "phase2":
        count = phase1(skel, enc, path)
    if mode == "phase1":
        json.dump({"skeleton": skel, "labelled_solutions": count}, sys.stdout, indent=1)
        return
    out = phase2(skel, enc, path)
    json.dump(out, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
