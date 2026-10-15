"""Plan C at n = 105 in the family F2MK: the router graph H is fixed per instance and only the
cross design (port assignment of the 16 units) and the pendant map (which doubles hang on which
port-like router) are free.

Instances: B = "K4,4 plus one hub-like vertex": remove two edges of K4,4 and join a ninth vertex to
their four endpoints (two disjoint edges -> simple B; two edges sharing a vertex -> a double edge
in B, realised by two port-like routers between the same hubs), then choose the 2k - 16 = 2 direct
(unsubdivided) edges among the 18 edges of B (no two parallel). H = the subdivision of the other
16 edges. Instances are deduplicated by the pynauty certificate of H with the colouring
{hub-like, port-like}; each is solved by twoblock_sat.Model(105, "F2MK") under assumptions fixing H.

Usage: python plan_c_105.py [CAP_SECONDS_PER_INSTANCE] [MAX_INSTANCES]
Output: runs/plan_c_105.jsonl (one record per instance), runs/plan_c_105_summary.json
"""
import collections
import itertools
import json
import pathlib
import sys
import time

import pynauty

import twoblock_sat as TS

HERE = pathlib.Path(__file__).resolve().parent
K = 9  # hub-like routers at n = 105


def k44_edges():
    return [(l, 4 + r) for l in range(4) for r in range(4)]  # hubs 0-3 = L, 4-7 = R


def base_multigraphs():
    """B as a list of (multi)edges on vertices 0..8 (8 = the new hub), for the two removal types."""
    E = k44_edges()
    out = {}
    for e1, e2 in itertools.combinations(E, 2):
        edges = [e for e in E if e not in (e1, e2)]
        for v in e1 + e2:
            edges.append((min(v, 8), max(v, 8)))
        shared = len(set(e1) & set(e2)) > 0
        out.setdefault(shared, []).append(sorted(edges))
    return out


def subdivided_H(edges, direct):
    """H on 9 hub-like (0..8) + 16 port-like (9..24) routers; `direct` = set of indices into edges."""
    H = set()
    nxt = K
    for i, (a, b) in enumerate(edges):
        if i in direct:
            H.add((a, b))
        else:
            H.add((a, nxt))
            H.add((b, nxt))
            nxt += 1
    assert nxt == K + 16
    return H


def certificate(H):
    adj = collections.defaultdict(set)
    for a, b in H:
        adj[a].add(b)
        adj[b].add(a)
    g = pynauty.Graph(K + 16, adjacency_dict={v: list(adj[v]) for v in adj},
                      vertex_coloring=[set(range(K)), set(range(K, K + 16))])
    return pynauty.certificate(g)


def instances():
    seen = {}
    for shared, blist in base_multigraphs().items():
        for edges in blist:
            for direct in itertools.combinations(range(len(edges)), 2):
                i, j = direct
                if edges[i] == edges[j]:
                    continue  # parallel direct edges would be a double edge
                H = subdivided_H(edges, set(direct))
                cert = certificate(H)
                if cert not in seen:
                    seen[cert] = {"shared_removal": shared, "B_edges": edges, "direct": [edges[i], edges[j]], "H": sorted(H)}
    return list(seen.values())


def main():
    cap = float(sys.argv[1]) if len(sys.argv) > 1 else 120
    maxinst = int(sys.argv[2]) if len(sys.argv) > 2 else 10 ** 9
    inst = instances()
    print(f"{len(inst)} non-isomorphic H instances", file=sys.stderr, flush=True)
    outp = HERE / "runs/plan_c_105.jsonl"
    done = set()
    if outp.exists():
        for line in outp.read_text().splitlines():
            done.add(json.loads(line)["index"])
    summary = collections.Counter()
    t_all = time.time()
    for idx, rec in enumerate(inst):
        if idx in done or idx >= maxinst:
            continue
        # router_sym=False: with H fixed, the model's router-relabelling rows are not WLOG, so only
        # "W[0] holds block-1 port (0,0)" (a relabelling of the doubles) is kept
        model = TS.Model(105, "F2MK", router_sym=False)
        R0 = model.W[32:]
        HUBR, PORTR = R0[:K], R0[K:]
        pos = {}
        for i in range(K):
            pos[i] = HUBR[i]
        for i in range(16):
            pos[K + i] = PORTR[i]
        Hset = set((min(pos[a], pos[b]), max(pos[a], pos[b])) for a, b in rec["H"])
        lits = []
        for u, v in itertools.combinations(R0, 2):
            var = model.evar(u, v)
            lits.append(var if (u, v) in Hset else -var)
        t0 = time.time()
        res = model.solve(cap=cap, assumptions=lits)
        out = {"index": idx, "shared_removal": rec["shared_removal"], "direct": rec["direct"],
               "result": {True: "SAT", False: "UNSAT", None: "TIMEOUT"}[res], "solve_time": time.time() - t0,
               "conflicts": getattr(model, "conflicts", None)}
        if res:
            edges = model.decode()
            out["bfs"] = TS.bfs_check(edges, 105)
            path = HERE / f"runs/105_4_4_planc_{idx}.edges"
            path.write_text("".join(f"{u} {v}\n" for u, v in edges))
            out["edges_file"] = str(path)
        with open(outp, "a") as f:
            f.write(json.dumps(out) + "\n")
        summary[out["result"]] += 1
        print(idx, out["result"], f"{out['solve_time']:.1f}s", file=sys.stderr, flush=True)
        del model
    json.dump({"instances": len(inst), "run": dict(summary), "wall": time.time() - t_all},
              open(HERE / "runs/plan_c_105_summary.json", "w"), indent=1)
    print(json.dumps({"instances": len(inst), "run": dict(summary)}))


if __name__ == "__main__":
    main()
