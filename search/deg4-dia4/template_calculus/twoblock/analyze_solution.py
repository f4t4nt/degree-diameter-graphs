"""Describe a decoded SAT solution: BFS check, number of subdivided K4,4 blocks, two-block statistics
of W (profiles, routers, B2, cross incidences), nauty certificate match against the 36 known 104s."""
import collections
import hashlib
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import verify  # noqa: E402
import known104_twoblock as K  # noqa: E402


def certificate(edges):
    import pynauty
    verts = sorted({v for e in edges for v in e})
    idx = {v: i for i, v in enumerate(verts)}
    adj = collections.defaultdict(list)
    for a, b in edges:
        adj[idx[a]].append(idx[b])
        adj[idx[b]].append(idx[a])
    g = pynauty.Graph(len(verts), adjacency_dict=dict(adj))
    return pynauty.certificate(g), pynauty.autgrp(g)[1]


def main():
    path = sys.argv[1]
    edges = verify.load_edges_file(path)
    verts, adj, problems = verify.build(edges)
    rec = {"file": path, "sha256": hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest(), "problems": problems}
    rec.update(verify.basic_checks(verts, adj))
    blocks = verify.k44_blocks(verts, adj)
    rec["k44_blocks"] = len(blocks)
    keys = list(blocks)
    rec["three_block_decompositions"] = len(verify.decompositions(blocks, len(verts)))
    views = []
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            if not (a & b):
                views.append(K.two_block_view(verts, adj, a, b, blocks, None))
    rec["two_block_views"] = views
    try:
        cert, aut = certificate(edges)
        rec["aut"] = aut
        known = {}
        for name, kedges in K.load_all():
            known[certificate(kedges)[0]] = name
        rec["isomorphic_known_class"] = known.get(cert)
    except ImportError:
        rec["isomorphic_known_class"] = "pynauty unavailable"
    print(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
