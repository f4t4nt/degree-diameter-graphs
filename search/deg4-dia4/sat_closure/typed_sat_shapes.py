"""Typed SAT driver for the shape lists of 106, 107 and 108 (data/shapes*.json[.gz]).

The model is the validated encoder of typed_slack_sat.py (imported unchanged: same variables,
clauses, redundant port-lemma and exactness clauses, symmetry breaking, tight mode for budget-768
shapes). The driver only adds:
  - a per-instance cap enforced with CaDiCaL conflict budgets in chunks (solve_limited is called
    repeatedly with conf_budget(CHUNK) until it decides or the wall-clock cap is exceeded -> TIMEOUT);
  - an order of expected difficulty: shapes with a gadget connector whose carried vector
    f != (2,2,2) first, unit-like shapes (every gadget connector has f = (2,2,2)) last;
  - results appended to the output JSON (and one record per line to OUT + "l") after every
    instance, so a batch can be resumed; a SAT answer is decoded, checked by independent BFS and
    written as runs/solution_<n>_idx<i>.edges;
  - --staged: hub-to-port / hub-to-hub model first, full model only if that is not UNSAT.

Usage: typed_sat_shapes.py SHAPES.json[.gz] --out=runs/NAME.json [--cap=SECONDS] [--idx=i,j,...]
                           [--order=easy|file] [--skip-unitlike] [--only-unitlike] [--chunk=CONFLICTS]
                           [--staged] [--extra=OTHER_RESULTS.json]
"""
import gzip
import hashlib
import json
import pathlib
import sys
import time

from pysat.solvers import Cadical153

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import typed_slack_sat as T  # noqa: E402


def unit_like(shape):
    return all(tuple(f) == (2, 2, 2) for g in shape["gadgets"] for f in g["f"])


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


def run_shape(shape, idx, cap, chunk, outdir, staged=False):
    """staged: stage 1 = hub-to-port / hub-to-hub model (full=False; UNSAT there is UNSAT for the
    graph problem), stage 2 = full model when stage 1 is SAT or undecided."""
    t0 = time.time()
    tight = shape.get("budget") == 768
    m = T.Model(shape["units"], shape["gadgets"], exact_units=True, lemma=True, full=not staged, tight=tight, exact_all=True)
    rec = {"idx": idx, "units": shape["units"], "potential": shape["potential"], "budget": shape.get("budget"),
           "gadgets": [(g["types"], g["edges"]) for g in shape["gadgets"]], "unit_like": unit_like(shape), "tight": tight,
           "vars": m.pool.top, "clauses": len(m.cl), "build_s": round(time.time() - t0, 1)}
    r, el, model, conflicts = solve_chunked(m.cl, cap, chunk)
    if staged:
        rec["stage1"] = {True: "SAT", False: "UNSAT", None: "TIMEOUT"}[r]
        rec["stage1_s"] = round(el, 1)
        if r is not False:
            t1 = time.time()
            m = T.Model(shape["units"], shape["gadgets"], exact_units=True, lemma=True, full=True, tight=tight, exact_all=True)
            rec["stage2_vars"], rec["stage2_clauses"], rec["stage2_build_s"] = m.pool.top, len(m.cl), round(time.time() - t1, 1)
            r, el, model, conflicts = solve_chunked(m.cl, cap, chunk)
    rec["result"] = {True: "SAT", False: "UNSAT", None: "TIMEOUT"}[r]
    rec["solve_s"] = round(el, 1)
    rec["conflict_budget_used"] = conflicts
    if r:
        edges = m.decode(model)
        rec["check"] = T.bfs_check(edges)
        n = rec["check"]["n"]
        path = outdir / f"solution_{n}_idx{idx}.edges"
        path.write_text("".join(f"{a} {b}\n" for a, b in edges))
        rec["edges_file"] = str(path)
        rec["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return rec


def main(args):
    opener = gzip.open if args[0].endswith(".gz") else open
    shapes = json.load(opener(args[0], "rt"))["shapes"]
    idxs = next((a.split("=")[1] for a in args if a.startswith("--idx=")), None)
    idxs = [int(x) for x in idxs.split(",")] if idxs else list(range(len(shapes)))
    cap = float(next((a.split("=")[1] for a in args if a.startswith("--cap=")), 1800))
    chunk = int(next((a.split("=")[1] for a in args if a.startswith("--chunk=")), 20000))
    out = next((a.split("=")[1] for a in args if a.startswith("--out=")), None)
    order = next((a.split("=")[1] for a in args if a.startswith("--order=")), "easy")
    if "--skip-unitlike" in args:
        idxs = [i for i in idxs if not unit_like(shapes[i])]
    if "--only-unitlike" in args:
        idxs = [i for i in idxs if unit_like(shapes[i])]
    if order == "easy":
        idxs.sort(key=lambda i: (unit_like(shapes[i]), -len(shapes[i]["gadgets"]), shapes[i]["units"], i))
    outdir = HERE / "runs"
    results = []
    # resume: records already in OUT (json list) or in OUT + "l" (one JSON record per line, written
    # after every instance; re-dumping the whole list each time cost 0.6 s per instance at 900 records)
    if out and pathlib.Path(out).exists():
        results = json.load(open(out))
    if out and pathlib.Path(out + "l").exists():
        seen_idx = {r["idx"] for r in results}
        for line in open(out + "l"):
            if line.strip():
                r = json.loads(line)
                if r["idx"] not in seen_idx:
                    results.append(r)
                    seen_idx.add(r["idx"])
    done = {r["idx"] for r in results}
    idxs = [i for i in idxs if i not in done]
    extra = [a for a in args if a.startswith("--extra=")]
    if extra:  # merge results of another file (e.g. a previous run's json) into the resume set
        for r in json.load(open(extra[0].split("=")[1])):
            if r["idx"] not in done:
                results.append(r)
                done.add(r["idx"])
        idxs = [i for i in idxs if i not in done]
    staged = "--staged" in args
    t_all = time.time()
    jl = open(out + "l", "a") if out else None
    for i in idxs:
        rec = run_shape(shapes[i], i, cap, chunk, outdir, staged=staged)
        rec["elapsed_total_s"] = round(time.time() - t_all, 1)
        print(json.dumps({k: rec.get(k) for k in ("idx", "units", "potential", "budget", "unit_like", "stage1", "stage1_s", "result", "solve_s", "vars", "clauses")}),
              file=sys.stderr, flush=True)
        if rec["result"] == "SAT":
            print("!!! SAT: " + json.dumps(rec), file=sys.stderr, flush=True)
        results.append(rec)
        if jl:
            jl.write(json.dumps(rec) + "\n")
            jl.flush()
    if out:
        json.dump(results, open(out, "w"), indent=1)
    else:
        json.dump(results, sys.stdout, indent=1)


if __name__ == "__main__":
    main(sys.argv[1:])
