"""Merge every typed-SAT verdict file for n = 108 against the complete shape lists (768:
shapes108_full_768.json, 11,163 shapes; 772: shapes108_full_772_noexcess.json, 7,056 shapes, of which
the 2,325 of shapes108_full_772_refined.json pass the refined excess matching).

This is the bookkeeping step of the 108 closure, kept as the record of how the per-shape verdicts of
several batches (local runs of typed_sat_shapes.py on the full lists and on the earlier S8 sub-family
lists, shards run on other machines, and the two all-triangle shapes decided by the skeleton method)
were combined into one verdict per shape. All records carry an "idx" into the shape file they were
run on; shapes are mapped between files by the canonical key (units, shapes_fast.shape_key(gadgets)).
The verdict files themselves are not part of this directory; the script expects them under runs/
(local batches) and under --shards=DIR/results/<batch>/typed.jsonl (shards).

Merged verdict per shape: SAT if any record says SAT (reported loudly), else UNSAT if any record says
UNSAT, else TIMEOUT if any record exists, else NOT_RUN.

Outputs (runs/): verdicts108_768.jsonl, verdicts108_772.jsonl (one line per shape of the full list:
idx, budget, units, unit_like, verdict, sources), merge_summary.json, open768_idx.txt, open772_refined_idx.txt,
open772_ineq_idx.txt (indices into the noexcess file of the inequality-only shapes without a verdict).

Usage: merge108e.py [--shards=DIR]
"""
import collections
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import shapes_fast as SF  # noqa: E402

SHARDS = HERE / "shards"
D = HERE / "runs"
C = HERE / "runs"


def load_shapes(path):
    return json.load(open(path))["shapes"]


def keys_of(shapes):
    return [(s["units"], SF.shape_key(s["gadgets"])) for s in shapes]


def unit_like(shape):
    return all(tuple(f) == (2, 2, 2) for g in shape["gadgets"] for f in g["f"])


def read_jsonl(path):
    return [json.loads(l) for l in open(path) if l.strip()]


def main(args):
    shards = pathlib.Path(next((a.split("=")[1] for a in args if a.startswith("--shards=")), SHARDS))
    out = HERE / "runs"
    out.mkdir(exist_ok=True)
    summary = {}

    full = {768: load_shapes(D / "shapes108_full_768.json"), 772: load_shapes(D / "shapes108_full_772_noexcess.json")}
    refined = load_shapes(D / "shapes108_full_772_refined.json")
    fkeys = {b: keys_of(full[b]) for b in full}
    rkeys = keys_of(refined)
    for b in full:
        assert all(s["budget"] == b for s in full[b]), b
        dup = len(fkeys[b]) - len(set(fkeys[b]))
        summary[f"full{b}"] = {"shapes": len(full[b]), "duplicate_keys": dup}
    assert all(s["budget"] == 772 for s in refined)
    index = {b: {k: i for i, k in enumerate(fkeys[b])} for b in full}
    # refined 772 -> noexcess 772
    ref2full = [index[772].get(k) for k in rkeys]
    summary["refined772"] = {"shapes": len(refined), "duplicate_keys": len(rkeys) - len(set(rkeys)),
                             "not_in_noexcess": sum(1 for i in ref2full if i is None)}
    assert summary["refined772"]["not_in_noexcess"] == 0, "refined list is not a subset of the noexcess list"
    refined_full_idx = set(ref2full)

    # S8 lists
    s8 = {"refined": load_shapes(C / "shapes108s8_refined.json"), "noexcess": load_shapes(C / "shapes108s8_noexcess.json"),
          "m8": load_shapes(C / "shapes108s8_m8_noexcess.json")}
    s8keys = {n: keys_of(v) for n, v in s8.items()}
    s8map = {}   # (list, idx) -> (budget, full idx)
    miss = collections.Counter()
    for n, v in s8.items():
        for i, s in enumerate(v):
            b = s["budget"]
            j = index[b].get(s8keys[n][i])
            if j is None:
                miss[n] += 1
            else:
                s8map[(n, i)] = (b, j)
    summary["s8_lists"] = {n: {"shapes": len(v), "not_in_full_list": miss[n]} for n, v in s8.items()}
    # S8 shapes in the full lists (by key, over all three S8 files)
    s8_full = {b: set() for b in full}
    for (n, i), (b, j) in s8map.items():
        s8_full[b].add(j)
    summary["s8_in_full"] = {b: len(v) for b, v in s8_full.items()}

    # verdict records: (budget, full idx) -> list of (source, result, solve_s, src idx)
    rec = collections.defaultdict(list)

    def add_full(path, b, mapper=None):
        n = 0
        for r in read_jsonl(path):
            assert r["budget"] == b, (path, r["idx"], r["budget"])
            j = r["idx"] if mapper is None else mapper[r["idx"]]
            rec[(b, j)].append((str(path), r["result"], r.get("solve_s"), r["idx"]))
            n += 1
        return n

    src_counts = {}
    src_counts[str(D / "typed108_full_768.jsonl")] = add_full(D / "typed108_full_768.jsonl", 768)
    for k in range(8):
        p = shards / "results" / f"sat108_768_s{k}of8" / "typed.jsonl"
        src_counts[str(p)] = add_full(p, 768)
    for p in sorted((shards / "results").glob("sat108e_768_i*/typed.jsonl")):   # task 2: one job per open shape
        src_counts[str(p)] = add_full(p, 768)
    for p in sorted((HERE / "runs").glob("typed108e_*.jsonl")):                # task 2: local runs of this directory
        src_counts[str(p)] = add_full(p, 768)
    src_counts[str(D / "typed108_full_772_refined.jsonl")] = add_full(D / "typed108_full_772_refined.jsonl", 772, ref2full)
    for k in range(4):
        p = shards / "results" / f"sat108_772r_s{k}of4" / "typed.jsonl"
        src_counts[str(p)] = add_full(p, 772, ref2full)

    def add_s8(path, listname):
        n = 0
        for r in read_jsonl(path):
            b, j = s8map[(listname, r["idx"])]
            assert r["budget"] == b
            rec[(b, j)].append((str(path), r["result"], r.get("solve_s"), r["idx"]))
            n += 1
        return n

    src_counts[str(C / "typed108s8.jsonl")] = add_s8(C / "typed108s8.jsonl", "refined")
    src_counts[str(C / "typed108s8_extra772.jsonl")] = add_s8(C / "typed108s8_extra772.jsonl", "noexcess")
    src_counts[str(C / "typed108s8_m8.jsonl")] = add_s8(C / "typed108s8_m8.jsonl", "m8")
    src_counts[str(D / "typed108s8_open.jsonl")] = add_s8(D / "typed108s8_open.jsonl", "m8")
    # skeleton method: the two all-triangle S8 shapes
    skel = json.load(open(D / "sat108_all.json"))  # output of `sat106.py 108 all`
    assert len(skel) == 54 and all(r["sat"] is False for r in skel), "sat108_all.json is not 54 x UNSAT"
    for i in (2011, 2015):
        s = s8["m8"][i]
        assert s["budget"] == 768 and s["units"] == 12 and len(s["gadgets"]) == 4 and all(len(g["types"]) == 3 and len(g["edges"]) == 3 for g in s["gadgets"]), i
        b, j = s8map[("m8", i)]
        rec[(b, j)].append((str(D / "sat108_all.json"), "UNSAT", None, i))
    summary["source_record_counts"] = src_counts

    # the 10 S8 shapes that were still open after the earlier batches
    open10 = [944, 996, 1033, 1330, 1333, 1334, 1335, 1380, 1381, 1795]
    summary["open_s8_10"] = {str(i): {"budget": s8map[("m8", i)][0], "full_idx": s8map[("m8", i)][1],
                                      "records": [(pathlib.Path(p).parent.parent.name + "/" + pathlib.Path(p).name, res, t) for (p, res, t, _) in rec[s8map[("m8", i)]]]}
                             for i in open10}

    # merged verdicts
    order = {"SAT": 0, "UNSAT": 1, "TIMEOUT": 2}
    per_budget = {}
    for b in full:
        rows = []
        cnt = collections.Counter()
        undecided = []
        for j, s in enumerate(full[b]):
            rs = rec.get((b, j), [])
            if not rs:
                v = "NOT_RUN"
            else:
                v = min((r[1] for r in rs), key=lambda x: order[x])
            cls = "S8" if j in s8_full[b] else ("refined" if (b == 772 and j in refined_full_idx) else ("ineq_only" if b == 772 else "new"))
            cnt[(cls, v)] += 1
            if v != "UNSAT":
                undecided.append(j)
            rows.append({"idx": j, "budget": b, "units": s["units"], "unit_like": unit_like(s), "class": cls, "verdict": v,
                         "sources": [{"file": str(pathlib.Path(p).relative_to(HERE) if str(p).startswith(str(HERE)) else pathlib.Path(p).relative_to(shards.parent)),
                                      "result": res, "solve_s": t, "src_idx": si} for (p, res, t, si) in rs]})
        with open(out / f"verdicts108_{b}.jsonl", "w") as f:
            for r in rows:
                f.write(json.dumps(r) + "\n")
        tot = collections.Counter(r["verdict"] for r in rows)
        per_budget[b] = {"shapes": len(rows), "verdicts": dict(tot), "by_class_verdict": {f"{c}/{v}": n for (c, v), n in sorted(cnt.items())},
                         "undecided_idx": undecided,
                         "sat_idx": [r["idx"] for r in rows if r["verdict"] == "SAT"],
                         "multi_source": sum(1 for r in rows if len(r["sources"]) > 1),
                         "conflicting_sources": [r["idx"] for r in rows if len({s["result"] for s in r["sources"]} & {"SAT", "UNSAT"}) > 1]}
        if b == 768:
            (out / "open768_idx.txt").write_text(",".join(map(str, undecided)) + "\n")
        else:
            ref_open = [j for j in undecided if j in refined_full_idx or j in s8_full[b]]
            ineq_open = [j for j in undecided if j not in refined_full_idx and j not in s8_full[b]]
            (out / "open772_refined_idx.txt").write_text(",".join(map(str, ref_open)) + "\n")
            (out / "open772_ineq_idx.txt").write_text(",".join(map(str, ineq_open)) + "\n")
            per_budget[b]["undecided_refined_or_s8"] = ref_open
            per_budget[b]["undecided_ineq_only_count"] = len(ineq_open)
            del per_budget[b]["undecided_idx"]
            # refined-list view: the 2,325 refined shapes (the shard indices are into the refined file)
            rv = collections.Counter()
            for i, j in enumerate(ref2full):
                rv[rows[j]["verdict"]] += 1
            per_budget[b]["refined_list_verdicts"] = dict(rv)
            per_budget[b]["refined_list_undecided_refined_idx"] = [i for i, j in enumerate(ref2full) if rows[j]["verdict"] != "UNSAT"]
    summary["per_budget"] = per_budget
    json.dump(summary, open(out / "merge_summary.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main(sys.argv[1:])
