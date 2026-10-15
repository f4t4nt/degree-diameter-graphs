"""Cube-and-conquer fallback for a typed-SAT shape that times out as one instance.

The model is the validated encoder typed_slack_sat.py (imported unchanged, same options
as typed_sat_shapes.py: exact_units, lemma, full, tight for budget 768, exact_all);
a cube adds unit clauses on the port variables A[(c, q)] of gadget connectors. The 16 cubes
"gadget connector C holds port (B, l, r)" for l, r in 0..3 (C with exactly one block-B port, so the
encoder's exactly-x_B(C) cardinality makes them a partition of the solution space) are a sound
case split: UNSAT on all 16 is UNSAT for the shape. A SAT answer is decoded, BFS-checked and written
next to --out like the other drivers.

Usage:
  cube_sat108e.py SHAPES.json --idx=I --cube=C:B,l,r[;C:B,l,r...] --cap=S --out=runs/NAME.json
  cube_sat108e.py SHAPES.json --idx=I --split=C:B --cap=S --out=runs/NAME.json     (all 16 cubes, one process)
  cube_sat108e.py SHAPES.json --idx=I --list                                        (gadget connectors and types)
"""
import hashlib
import json
import pathlib
import sys
import time

from pysat.solvers import Cadical153

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import typed_slack_sat as T  # noqa: E402


def solve_chunked(clauses, cap, chunk):
    s = Cadical153(bootstrap_with=clauses)
    t0 = time.time()
    r = None
    conflicts = 0
    while True:
        s.conf_budget(chunk)
        r = s.solve_limited()
        conflicts += chunk
        if r is not None or time.time() - t0 > cap:
            break
    el = time.time() - t0
    model = s.get_model() if r else None
    s.delete()
    return r, el, model, conflicts


def parse_cube(spec):
    lits = []
    for part in spec.split(";"):
        if not part:
            continue
        neg = part.startswith("-")
        c, rest = part.lstrip("-").split(":")
        b, l, r = (int(x) for x in rest.split(","))
        lits.append((int(c), T.PIDX[(b, l, r)], not neg))
    return lits


def run_cube(shape, idx, cube, cap, chunk, outdir):
    t0 = time.time()
    tight = shape.get("budget") == 768
    m = T.Model(shape["units"], shape["gadgets"], exact_units=True, lemma=True, full=True, tight=tight, exact_all=True)
    for (c, q, pos) in cube:
        assert m.types[c] is not None, f"connector {c} is a unit member"
        lit = m.A.get((c, q))
        if lit is None:   # the connector has no port in that block: the cube is trivially UNSAT / trivially true
            if pos:
                m.cl.append([])  # empty clause: UNSAT
            continue
        m.cl.append([lit if pos else -lit])
    rec = {"idx": idx, "units": shape["units"], "potential": shape["potential"], "budget": shape.get("budget"),
           "cube": [(c, T.PORTS[q], pos) for (c, q, pos) in cube], "tight": tight, "vars": m.pool.top, "clauses": len(m.cl),
           "build_s": round(time.time() - t0, 1)}
    if any(len(cl) == 0 for cl in m.cl):
        rec.update({"result": "UNSAT", "solve_s": 0.0, "conflict_budget_used": 0, "note": "connector has no port in the cube's block"})
        return rec
    r, el, model, conflicts = solve_chunked(m.cl, cap, chunk)
    rec["result"] = {True: "SAT", False: "UNSAT", None: "TIMEOUT"}[r]
    rec["solve_s"] = round(el, 1)
    rec["conflict_budget_used"] = conflicts
    if r:
        edges = m.decode(model)
        rec["check"] = T.bfs_check(edges)
        n = rec["check"]["n"]
        path = outdir / f"solution_{n}_idx{idx}_cube.edges"
        path.write_text("".join(f"{a} {b}\n" for a, b in edges))
        rec["edges_file"] = str(path)
        rec["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return rec


def main(args):
    shapes = json.load(open(args[0]))["shapes"]
    idx = int(next(a.split("=")[1] for a in args if a.startswith("--idx=")))
    shape = shapes[idx]
    if "--list" in args:
        k = 0
        for gi, g in enumerate(shape["gadgets"]):
            for t, f in zip(g["types"], g["f"]):
                print(k, "gadget", gi, "type", t, "f", f)
                k += 1
        return
    cap = float(next((a.split("=")[1] for a in args if a.startswith("--cap=")), 3300))
    chunk = int(next((a.split("=")[1] for a in args if a.startswith("--chunk=")), 20000))
    out = next(a.split("=")[1] for a in args if a.startswith("--out="))
    outdir = pathlib.Path(out).resolve().parent
    cube = next((a.split("=")[1] for a in args if a.startswith("--cube=")), None)
    split = next((a.split("=")[1] for a in args if a.startswith("--split=")), None)
    cubes = []
    if cube is not None:
        cubes.append(parse_cube(cube))
    if split is not None:
        c, b = (int(x) for x in split.split(":"))
        cubes += [[(c, T.PIDX[(b, l, r)], True)] for l in range(4) for r in range(4)]
    t_all = time.time()
    with open(out + "l", "a") as jl:
        for cb in cubes:
            rec = run_cube(shape, idx, cb, cap, chunk, outdir)
            rec["elapsed_total_s"] = round(time.time() - t_all, 1)
            print(json.dumps({k: rec.get(k) for k in ("idx", "cube", "result", "solve_s", "vars", "clauses")}), file=sys.stderr, flush=True)
            if rec["result"] == "SAT":
                print("!!! SAT: " + json.dumps(rec), file=sys.stderr, flush=True)
            jl.write(json.dumps(rec) + "\n")
            jl.flush()


if __name__ == "__main__":
    main(sys.argv[1:])
