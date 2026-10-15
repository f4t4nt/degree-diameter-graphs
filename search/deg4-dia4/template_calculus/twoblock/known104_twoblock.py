"""Calibration on the 36 known 104s: how many subdivided K4,4 blocks does each contain, and what does
each graph look like in the two-block view (two of its blocks kept, the rest of the graph = W)?

For every graph: all subdivided K4,4 blocks (verify.k44_blocks, exhaustive), the number of
three-block decompositions, and for each of the three block pairs the two-block statistics:
profile histogram of W, routers, cross incidences sum x1 x2, the two-block budget
B2 = 2 sum_c (x1 f2 + x2 f1), per-port f-sums, per-hub foreign-port reach and the W-W / W-port /
hub-W conditions in the G'' formulation (checked by BFS, so this also validates the pair conditions
of RESULTS.md section 1 against the real distances).

Output: runs/known104_twoblock.json
"""
import collections
import itertools
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import verify  # noqa: E402


def load_all():
    corpus = verify.load_corpus()
    out = []
    for i, rep in enumerate(corpus["representatives"]):
        out.append((f"class{i}", [tuple(e) for e in rep["edges_sorted"]]))
    return out


def two_block_view(verts, adj, b1, b2, blocks, dist):
    """Statistics of the graph seen as blocks b1, b2 plus remainder W."""
    hub_block, port_block, port_hubs = {}, {}, {}
    for i, b in enumerate((b1, b2)):
        L, R, pmap = blocks[b]
        for v in L + R:
            hub_block[v] = i
        for (l, r), q in pmap.items():
            port_block[q] = i
            port_hubs[q] = (l, r)
    W = [v for v in verts if v not in hub_block and v not in port_block]
    Wset = set(W)
    x = {c: [0, 0] for c in W}
    plain = {c: set() for c in W}
    for c in W:
        for u in adj[c]:
            if u in port_block:
                x[c][port_block[u]] += 1
            elif u in Wset:
                plain[c].add(u)
            else:
                raise AssertionError("W vertex adjacent to a hub")
    f = {c: [x[c][j] + sum(x[d][j] for d in plain[c]) for j in range(2)] for c in W}
    T = {c: x[c][0] * f[c][1] + x[c][1] * f[c][0] for c in W}
    B2 = 2 * sum(T.values())
    cross = sum(x[c][0] * x[c][1] for c in W)
    hist = collections.Counter(tuple(x[c]) for c in W)
    routers = [c for c in W if sum(x[c]) == 0]
    # per-port f-sum of the holders in the foreign coordinate
    port_fsum = {}
    for q, j in port_block.items():
        holders = [u for u in adj[q] if u in Wset]
        assert len(holders) == 2
        port_fsum[q] = sum(f[h][1 - j] for h in holders)
    # per-hub foreign port reach S_j(a) with multiplicity and distinct count
    hub_reach = {}
    for a, i in hub_block.items():
        j = 1 - i
        multiset = collections.Counter()
        for p in adj[a]:
            for c in adj[p]:
                if c in Wset:
                    for q in adj[c]:
                        if q in port_block and port_block[q] == j:
                            multiset[q] += 1
                    for d in plain[c]:
                        for q in adj[d]:
                            if q in port_block and port_block[q] == j:
                                multiset[q] += 1
        hub_reach[a] = (sum(multiset.values()), len(multiset))
    # routers: port incidences within plain radius 2 (hub condition) and pairwise distances
    router_ball_ports = []
    for r in routers:
        ball = set(plain[r])
        for d in list(ball):
            ball |= plain[d]
        ball.discard(r)
        hubs = set()
        for c in ball:
            for q in adj[c]:
                if q in port_block:
                    hubs |= set(port_hubs[q])
        router_ball_ports.append((sum(sum(x[c]) for c in ball), len(hubs)))
    # routers' induced structure: is it the third block?
    rset = set(routers)
    router_deg = collections.Counter(len(plain[r] & rset) for r in routers)
    # star exemption check: W-W pairs at distance > 4 in G' (no hubs) must be star-mates
    return {"W": len(W), "profiles": {str(k): v for k, v in sorted(hist.items())},
            "routers": len(routers), "router_internal_degrees": dict(router_deg),
            "cross_incidences": cross, "B2": B2,
            "port_fsum_min": min(port_fsum.values()), "port_fsum_max": max(port_fsum.values()),
            "hub_reach_multiplicity_min": min(v[0] for v in hub_reach.values()),
            "hub_reach_multiplicity_max": max(v[0] for v in hub_reach.values()),
            "hub_reach_distinct_min": min(v[1] for v in hub_reach.values()),
            "router_ball_port_incidences_min": min(v[0] for v in router_ball_ports) if routers else None,
            "router_ball_hubs_covered_min": min(v[1] for v in router_ball_ports) if routers else None}


def main():
    results = []
    for name, edges in load_all():
        verts, adj, problems = verify.build(edges)
        assert not problems
        blocks = verify.k44_blocks(verts, adj)
        decs = verify.decompositions(blocks, len(verts))
        keys = list(blocks)
        dist = None
        rec = {"name": name, "n": len(verts), "k44_blocks": len(blocks), "three_block_decompositions": len(decs),
               "pairwise_disjoint_block_pairs": sum(1 for a, b in itertools.combinations(keys, 2) if not (a & b)),
               "views": []}
        for a, b in itertools.combinations(keys, 2):
            if a & b:
                continue
            rec["views"].append(two_block_view(verts, adj, a, b, blocks, dist))
        results.append(rec)
        v0 = rec["views"][0]
        print(name, "blocks", len(blocks), "decomps", len(decs), "views", len(rec["views"]),
              "B2", [v["B2"] for v in rec["views"]], "cross", [v["cross_incidences"] for v in rec["views"]],
              "routers", [v["routers"] for v in rec["views"]], v0["profiles"], file=sys.stderr, flush=True)
    summary = {"graphs": len(results),
               "k44_block_counts": dict(collections.Counter(r["k44_blocks"] for r in results)),
               "decomposition_counts": dict(collections.Counter(r["three_block_decompositions"] for r in results)),
               "B2_values": dict(collections.Counter(v["B2"] for r in results for v in r["views"])),
               "cross_incidence_values": dict(collections.Counter(v["cross_incidences"] for r in results for v in r["views"])),
               "router_counts": dict(collections.Counter(v["routers"] for r in results for v in r["views"])),
               "port_fsum_min": min(v["port_fsum_min"] for r in results for v in r["views"]),
               "port_fsum_max": max(v["port_fsum_max"] for r in results for v in r["views"]),
               "hub_reach_multiplicity": sorted({(v["hub_reach_multiplicity_min"], v["hub_reach_multiplicity_max"]) for r in results for v in r["views"]}),
               "profile_sets": sorted({tuple(sorted(v["profiles"])) for r in results for v in r["views"]}),
               "router_ball_port_incidences_min": min(v["router_ball_port_incidences_min"] for r in results for v in r["views"]),
               "router_ball_hubs_covered_min": min(v["router_ball_hubs_covered_min"] for r in results for v in r["views"])}
    json.dump({"summary": summary, "graphs": results}, open(HERE / "runs/known104_twoblock.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
