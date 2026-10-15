"""Simple connected cubic base graphs on K vertices for lifts over groups of
order M with target diameter D, prefiltered by the non-backtracking walk
condition and written in lift.c base format.

A lift of a simple base graph over a group of order M has diameter <= D only if
for every ordered pair (u, v) of base vertices the number of non-backtracking
walks of length <= D from u to v is >= M (every element of the fibre over v has
to be the net voltage of such a walk; u == v needs >= M - 1 closed walks of
length 1..D).  The script samples random cubic graphs (networkx), dedupes by
pynauty certificate, keeps those passing the condition (and, for --nonbip, the
non-bipartite ones) and writes them sorted by their walk-count slack.

usage: cubic_bases.py K M D SAMPLES OUTDIR [--nonbip] [--seed S]
"""
import argparse
import os
import sys

import networkx as nx
import numpy as np
import pynauty


def nb_walk_counts(G, D):
    """counts[u][v] = number of non-backtracking walks of length 1..D from u to v."""
    nodes = sorted(G.nodes)
    idx = {v: i for i, v in enumerate(nodes)}
    k = len(nodes)
    darts = []
    for u, v in G.edges:
        darts.append((idx[u], idx[v]))
        darts.append((idx[v], idx[u]))
    nd = len(darts)
    rev = {}
    for i, (a, b) in enumerate(darts):
        rev[i] = darts.index((b, a))
    out_of = [[i for i, (a, b) in enumerate(darts) if a == u] for u in range(k)]
    counts = np.zeros((k, k), dtype=np.int64)
    for s in range(k):
        cnt = np.zeros(nd, dtype=np.int64)
        for i in out_of[s]:
            cnt[i] = 1
        for ell in range(1, D + 1):
            for i, (a, b) in enumerate(darts):
                counts[s, b] += cnt[i]
            if ell == D:
                break
            new = np.zeros(nd, dtype=np.int64)
            for i, (a, b) in enumerate(darts):
                if cnt[i] == 0:
                    continue
                for j in out_of[b]:
                    if j != rev[i]:
                        new[j] += cnt[i]
            cnt = new
    return counts


def girth(G):
    best = None
    for c in nx.minimum_cycle_basis(G):
        if best is None or len(c) < best:
            best = len(c)
    return best


def write_base(path, G):
    nodes = sorted(G.nodes)
    idx = {v: i for i, v in enumerate(nodes)}
    k = len(nodes)
    T = nx.bfs_tree(G, nodes[0])
    tree = set(frozenset((idx[a], idx[b])) for a, b in T.edges)
    edges = [(min(idx[a], idx[b]), max(idx[a], idx[b])) for a, b in G.edges]
    slots = [e for e in edges if frozenset(e) in tree] + [e for e in edges if frozenset(e) not in tree]
    with open(path, "w") as f:
        f.write(f"{k} {len(slots)} {len(tree)}\n")
        for u, v in slots:
            f.write(f"{u} {v} 0\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("K", type=int)
    ap.add_argument("M", type=int)
    ap.add_argument("D", type=int)
    ap.add_argument("SAMPLES", type=int)
    ap.add_argument("OUTDIR")
    ap.add_argument("--nonbip", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--prefix", default=None)
    args = ap.parse_args()
    K, M, D = args.K, args.M, args.D
    os.makedirs(args.OUTDIR, exist_ok=True)
    seen = {}
    rejected_bip = rejected_walk = 0
    for t in range(args.SAMPLES):
        G = nx.random_regular_graph(3, K, seed=args.seed * 1000003 + t)
        if not nx.is_connected(G):
            continue
        g = pynauty.Graph(K, adjacency_dict={v: list(G.neighbors(v)) for v in G.nodes})
        cert = pynauty.certificate(g)
        if cert in seen:
            continue
        bip = nx.is_bipartite(G)
        if args.nonbip and bip:
            seen[cert] = None
            rejected_bip += 1
            continue
        counts = nb_walk_counts(G, D)
        off = counts + np.eye(K, dtype=np.int64)  # trivial walk for u == v
        slack = int(off.min()) - M
        if slack < 0:
            seen[cert] = None
            rejected_walk += 1
            continue
        _, a, b, _, _ = pynauty.autgrp(g)
        seen[cert] = dict(G=G, slack=slack, mean=float(off.mean()), girth=girth(G), diam=nx.diameter(G),
                          aut=int(round(a * 10 ** b)), bip=bip)
    kept = [v for v in seen.values() if v is not None]
    kept.sort(key=lambda r: (-r["slack"], -r["mean"]))
    prefix = args.prefix or f"c{K}"
    with open(os.path.join(args.OUTDIR, "index.txt"), "a") as idx:
        idx.write(f"# K={K} M={M} D={D} samples={args.SAMPLES} distinct={len(seen)} kept={len(kept)} "
                  f"rejected_bipartite={rejected_bip} rejected_walkcount={rejected_walk}\n")
        for i, r in enumerate(kept):
            name = f"{prefix}_{i:04d}"
            write_base(os.path.join(args.OUTDIR, name + ".txt"), r["G"])
            idx.write(f"{name} k={K} slack={r['slack']} meanwalks={r['mean']:.1f} girth={r['girth']} "
                      f"diam={r['diam']} aut={r['aut']} bipartite={r['bip']}\n")
    print(f"K={K} M={M} D={D}: {len(seen)} distinct connected cubic graphs sampled, kept {len(kept)} "
          f"(rejected: bipartite {rejected_bip}, walk count {rejected_walk})", flush=True)


if __name__ == "__main__":
    main()
