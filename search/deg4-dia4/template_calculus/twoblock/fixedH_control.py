"""F2MK model with the router graph H fixed by assumptions (H = subdivided K4,4 at n = 104): how
much faster is the free-remainder SAT when the router graph is given? Writes runs/fixedH_104.json
and, if SAT, runs/fixedH_104.edges (BFS-checked)."""
import json
import pathlib
import sys
import time

import twoblock_sat as TS

HERE = pathlib.Path(__file__).resolve().parent


def k44_router_graph(k=8):
    """Subdivided K4,4 on hub-like routers 0..7 (order L0, R0, L1, L2, L3, R1, R2, R3, so that port
    (0,0) is adjacent to hub-like 0 and 1) and port-like routers 8..23 (port (l, r) -> 8 + 4 l + r)."""
    hub_index = {("L", 0): 0, ("R", 0): 1, ("L", 1): 2, ("L", 2): 3, ("L", 3): 4, ("R", 1): 5, ("R", 2): 6, ("R", 3): 7}
    edges = set()
    for l in range(4):
        for r in range(4):
            p = 8 + 4 * l + r
            edges.add((hub_index[("L", l)], p))
            edges.add((hub_index[("R", r)], p))
    return edges


def main():
    n = 104
    cap = float(sys.argv[1]) if len(sys.argv) > 1 else 300
    model = TS.Model(n, "F2MK", router_sym=False)
    R0 = model.W[32:]
    H = k44_router_graph()
    lits = []
    for i in range(len(R0)):
        for j in range(i + 1, len(R0)):
            var = model.evar(R0[i], R0[j])
            lits.append(var if (i, j) in H else -var)
    t0 = time.time()
    res = model.solve(cap=cap, assumptions=lits)
    rec = {"n": n, "family": "F2MK", "H": "subdivided K4,4 (fixed by assumptions)", "cap": cap,
           "result": {True: "SAT", False: "UNSAT", None: "TIMEOUT"}[res], "solve_time": time.time() - t0,
           "conflict_budget_used": getattr(model, "conflicts", None)}
    if res:
        edges = model.decode()
        rec["bfs"] = TS.bfs_check(edges, n)
        path = HERE / "runs/fixedH_104.edges"
        path.write_text("".join(f"{u} {v}\n" for u, v in edges))
        rec["edges_file"] = str(path)
    json.dump(rec, open(HERE / "runs/fixedH_104.json", "w"), indent=1)
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
