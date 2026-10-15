"""Pair-to-triangle surgery on the known perfect 104s.

A matched pair of connectors {c, c'} owns two ports in each block (verified for all 35 classes).
Replacing the pair by three doubles d12 = (p, q), d23 = (q', r), d13 = (p', r') joined in a triangle
preserves every hub-to-port reach set (the "virtual pair") and adds one vertex. Converting t pairs
gives n = 104 + t; converting all 16 gives the 120-vertex shape the counting relaxation singles out.

Usage:
  python surgery.py t1            # every class, every pair, all 8 orientations (n = 105)
  python surgery.py t16 [LIMIT]   # every class: orientations with exact hub-to-hub coverage, n = 120
"""
import collections
import itertools
import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import verify  # noqa: E402


def load_classes():
    corpus = verify.load_corpus()
    out = []
    for i, rep in enumerate(corpus["representatives"]):
        edges = [tuple(e) for e in rep["edges_sorted"]]
        verts, adj, _ = verify.build(edges)
        blocks = verify.k44_blocks(verts, adj)
        trip = verify.decompositions(blocks, len(verts))[0]
        out.append((f"class{i}", verts, adj, trip, blocks))
    return out


class Template:
    def __init__(self, verts, adj, trip, blocks):
        self.adj = {v: set(adj[v]) for v in verts}
        self.n = len(verts)
        self.port_block, self.port_hubs, self.hub_block = {}, {}, {}
        for i, b in enumerate(trip):
            L, R, pmap = blocks[b]
            for v in L + R:
                self.hub_block[v] = i
            for (l, r), q in pmap.items():
                self.port_block[q] = i
                self.port_hubs[q] = (l, r)
        self.connectors = [v for v in verts if v not in self.port_block and v not in self.hub_block]
        cset = set(self.connectors)
        self.partner = {}
        self.ports_of = {}
        for c in self.connectors:
            self.ports_of[c] = [u for u in self.adj[c] if u in self.port_block]
            nb = [u for u in self.adj[c] if u in cset]
            assert len(nb) == 1 and len(self.ports_of[c]) == 3
            self.partner[c] = nb[0]
        self.pairs = sorted({tuple(sorted((c, self.partner[c]))) for c in self.connectors})
        self.pair_ports = {}
        for c, c2 in self.pairs:
            per = [[], [], []]
            for q in self.ports_of[c] + self.ports_of[c2]:
                per[self.port_block[q]].append(q)
            assert all(len(x) == 2 for x in per), per
            self.pair_ports[(c, c2)] = per

    def kind(self, v, new):
        if v in new:
            return "double"
        if v in self.hub_block:
            return "hub"
        if v in self.port_block:
            return "port"
        return "connector"

    def surgery(self, conversions):
        """conversions: list of (pair, orientation) with orientation bits (sp, sq, sr)."""
        adj = {v: set(s) for v, s in self.adj.items()}
        new = set()
        nxt = max(adj) + 1
        for (c, c2), (sp, sq, sr) in conversions:
            P1, P2, P3 = self.pair_ports[(c, c2)]
            p, p2 = (P1[0], P1[1]) if sp == 0 else (P1[1], P1[0])
            q, q2 = (P2[0], P2[1]) if sq == 0 else (P2[1], P2[0])
            r, r2 = (P3[0], P3[1]) if sr == 0 else (P3[1], P3[0])
            for x in (c, c2):
                for u in adj[x]:
                    adj[u].discard(x)
                del adj[x]
            d12, d23, d13 = nxt, nxt + 1, nxt + 2
            nxt += 3
            for d, (u, w) in ((d12, (p, q)), (d23, (q2, r)), (d13, (p2, r2))):
                adj[d] = {u, w}
                adj[u].add(d)
                adj[w].add(d)
            adj[d12] |= {d23, d13}
            adj[d23] |= {d12, d13}
            adj[d13] |= {d12, d23}
            new |= {d12, d23, d13}
        return adj, new

    def hh_rows(self, pair, orientation):
        """Hub-pair sets covered by the three doubles of an oriented triangle, per block pair."""
        P1, P2, P3 = self.pair_ports[pair]
        sp, sq, sr = orientation
        p, p2 = (P1[0], P1[1]) if sp == 0 else (P1[1], P1[0])
        q, q2 = (P2[0], P2[1]) if sq == 0 else (P2[1], P2[0])
        r, r2 = (P3[0], P3[1]) if sr == 0 else (P3[1], P3[0])
        rows = {}
        for key, (u, w) in (((0, 1), (p, q)), ((1, 2), (q2, r)), ((0, 2), (p2, r2))):
            rows[key] = {(a, b) for a in self.port_hubs[u] for b in self.port_hubs[w]}
        return rows


def defect(adj):
    verts = sorted(adj)
    idx = {v: i for i, v in enumerate(verts)}
    n = len(verts)
    A = np.eye(n, dtype=np.int64)
    for v in verts:
        for u in adj[v]:
            A[idx[v], idx[u]] = 1
    R2 = (A @ A) > 0
    R4 = (R2.astype(np.int64) @ R2.astype(np.int64)) > 0
    bad = [(verts[i], verts[j]) for i, j in zip(*np.where(~R4)) if i < j]
    degs = collections.Counter(len(adj[v]) for v in verts)
    return bad, degs


def run_t1(classes):
    results = []
    for name, verts, adj, trip, blocks in classes:
        T = Template(verts, adj, trip, blocks)
        best = None
        kinds = collections.Counter()
        hist = collections.Counter()
        for pair in T.pairs:
            for ori in itertools.product((0, 1), repeat=3):
                adj2, new = T.surgery([(pair, ori)])
                bad, degs = defect(adj2)
                assert set(degs) == {4} and len(adj2) == 105
                hist[len(bad)] += 1
                for u, v in bad:
                    kinds[tuple(sorted((T.kind(u, new), T.kind(v, new))))] += 1
                if best is None or len(bad) < best[0]:
                    best = (len(bad), pair, ori, bad)
        results.append({"class": name, "pairs": len(T.pairs), "min_defect": best[0],
                        "defect_histogram": dict(sorted(hist.items())),
                        "bad_pair_kinds_total": {f"{a}-{b}": v for (a, b), v in sorted(kinds.items())},
                        "best_example": {"pair": best[1], "orientation": best[2], "bad": best[3]}})
        print(name, "min defect", best[0], dict(sorted(hist.items())), file=sys.stderr)
    return results


def run_t16(classes, limit):
    results = []
    for name, verts, adj, trip, blocks in classes:
        T = Template(verts, adj, trip, blocks)
        rows = {pair: {ori: T.hh_rows(pair, ori) for ori in itertools.product((0, 1), repeat=3)}
                for pair in T.pairs}
        solutions = []
        nodes = [0]

        def dfs(i, cov, chosen):
            nodes[0] += 1
            if len(solutions) >= limit:
                return
            if i == len(T.pairs):
                solutions.append(list(chosen))
                return
            pair = T.pairs[i]
            for ori, rws in rows[pair].items():
                if all(not (rws[k] & cov[k]) for k in rws):
                    cov2 = {k: cov[k] | rws[k] for k in cov}
                    chosen.append((pair, ori))
                    dfs(i + 1, cov2, chosen)
                    chosen.pop()
                    if len(solutions) >= limit:
                        return

        dfs(0, {(0, 1): set(), (1, 2): set(), (0, 2): set()}, [])
        rec = {"class": name, "search_nodes": nodes[0], "exact_hh_orientations_found": len(solutions),
               "limit": limit}
        if solutions:
            best = None
            hist = collections.Counter()
            kinds = collections.Counter()
            for sol in solutions:
                adj2, new = T.surgery(sol)
                bad, degs = defect(adj2)
                assert set(degs) == {4} and len(adj2) == 120
                hist[len(bad)] += 1
                for u, v in bad:
                    kinds[tuple(sorted((T.kind(u, new), T.kind(v, new))))] += 1
                if best is None or len(bad) < best[0]:
                    best = (len(bad), sol, bad)
            rec.update({"min_defect_120": best[0], "defect_histogram": dict(sorted(hist.items())),
                        "bad_pair_kinds_total": {f"{a}-{b}": v for (a, b), v in sorted(kinds.items())},
                        "best_solution": best[1], "best_bad_pairs": best[2][:40]})
        print(name, rec.get("exact_hh_orientations_found"), "solutions,", "min defect",
              rec.get("min_defect_120"), file=sys.stderr)
        results.append(rec)
    return results


def main():
    mode = sys.argv[1]
    classes = load_classes()
    if mode == "t1":
        out = run_t1(classes)
    elif mode == "t16":
        out = run_t16(classes, int(sys.argv[2]) if len(sys.argv) > 2 else 2000)
    else:
        raise SystemExit("mode t1 or t16")
    print(json.dumps(out, indent=1, default=list))


if __name__ == "__main__":
    main()
