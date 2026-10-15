"""Check the exported class files: plain BFS (order, degrees, diameter, pairs beyond 4) and nauty
certificates (pairwise non-isomorphic). Lines starting with '#' are comments.

Usage: python check_classes.py class*.edges.txt
       python check_classes.py --match FILE... -- class*.edges.txt
           (also report, for every FILE, which class it is isomorphic to, if any)
"""
import collections
import sys

import pynauty


def load(path):
    edges = []
    for line in open(path):
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) == 2:
            edges.append((int(parts[0]), int(parts[1])))
    return edges


def bfs_stats(edges):
    adj = collections.defaultdict(set)
    for u, v in edges:
        assert u != v, "loop"
        assert v not in adj[u], "duplicate edge"
        adj[u].add(v)
        adj[v].add(u)
    verts = sorted(adj)
    assert verts == list(range(len(verts))), "vertices are not 0..n-1"
    n = len(verts)
    degs = sorted({len(adj[v]) for v in verts})
    diam = 0
    far = 0
    for s in verts:
        dist = {s: 0}
        frontier = [s]
        while frontier:
            nxt = []
            for u in frontier:
                for w in adj[u]:
                    if w not in dist:
                        dist[w] = dist[u] + 1
                        nxt.append(w)
            frontier = nxt
        assert len(dist) == n, "disconnected"
        diam = max(diam, max(dist.values()))
        far += sum(1 for t, d in dist.items() if t > s and d > 4)
    return n, degs, diam, far


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
    args = sys.argv[1:]
    match = []
    if args and args[0] == "--match":
        sep = args.index("--")
        match, args = args[1:sep], args[sep + 1:]
    certs = {}
    for path in args:
        edges = load(path)
        n, degs, diam, far = bfs_stats(edges)
        cert, aut = certificate(edges)
        dup = certs.get(cert)
        certs.setdefault(cert, path)
        print(f"{path}: n {n} edges {len(edges)} degrees {degs} diameter {diam} pairs beyond 4: {far} "
              f"automorphisms {int(aut)}" + (f"  ISOMORPHIC TO {dup}" if dup else ""))
    print(f"{len(args)} files, {len(certs)} distinct certificates")
    for path in match:
        cert, _ = certificate(load(path))
        print(f"{path}: isomorphic to {certs.get(cert, 'none of them')}")


if __name__ == "__main__":
    main()
