"""Reduce the labelled diameter-4 realisations of a skeleton (runs/all104_skel<k>_solutions.txt,
written by enumerate104.py) to isomorphism classes with nauty certificates (pynauty), and compare
them with the classes in ../classes/. Usage: python classify104.py <skeleton index>"""
import collections
import hashlib
import json
import pathlib
import sys

import pynauty

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "template_calculus"))
import sat105  # noqa: E402
import enumerate104  # noqa: E402
import verify  # noqa: E402


def certificate(edges):
    verts = sorted({v for e in edges for v in e})
    idx = {v: i for i, v in enumerate(verts)}
    adj = collections.defaultdict(list)
    for a, b in edges:
        adj[idx[a]].append(idx[b])
        adj[idx[b]].append(idx[a])
    g = pynauty.Graph(len(verts), adjacency_dict=dict(adj))
    return pynauty.certificate(g), pynauty.autgrp(g)[1]


def main():
    skel = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    units = sat105.load_skeletons()[skel]
    enc = sat105.Encoder(units, "104")
    classes = {}
    n = 0
    for line in open(f"runs/all104_skel{skel}_solutions.txt"):
        opts = [int(x) for x in line.split()]
        edges = enumerate104.edges_from_options(enc, opts)
        cert, aut = certificate(edges)
        if cert not in classes:
            v = sat105.verify_edges(edges)
            assert v["diameter"] == 4 and v["bad_pairs"] == 0 and v["degrees"] == {4: 104}, v
            classes[cert] = {"count": 0, "aut": aut, "edges": edges, "bent": None}
        classes[cert]["count"] += 1
        n += 1
    print("labelled:", n, "classes:", len(classes), file=sys.stderr, flush=True)
    corpus = verify.load_corpus()
    known_certs = {}
    for i, rep in enumerate(corpus["representatives"]):
        cert, aut = certificate([tuple(e) for e in rep["edges_sorted"]])
        known_certs[cert] = i
    recovered = [known_certs[c] for c in classes if c in known_certs]
    new = [c for c in classes if c not in known_certs]
    nk = len(corpus["representatives"])
    print("known recovered:", len(recovered), f"of {nk}; new classes:", len(new), file=sys.stderr, flush=True)
    out = {"labelled_solutions": n, "isomorphism_classes": len(classes), "known_recovered": sorted(recovered),
           "known_missing": sorted(set(range(nk)) - set(recovered)), "new_classes": len(new),
           "automorphism_group_orders": collections.Counter(int(c["aut"]) for c in classes.values()),
           "class_counts": sorted((c["count"] for c in classes.values()), reverse=True)}
    out["skeleton"] = skel
    with open(f"runs/all104_skel{skel}_classes.json", "w") as f:
        json.dump([{"edges": c["edges"], "aut": c["aut"], "labelled_count": c["count"], "known_index": known_certs.get(cert)}
                   for cert, c in classes.items()], f)
    json.dump(out, sys.stdout, indent=1, default=str)


if __name__ == "__main__":
    main()
