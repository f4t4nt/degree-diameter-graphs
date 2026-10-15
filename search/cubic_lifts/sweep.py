"""Sequential sweep driver for lift: for each target order n, each k | n with
base graphs available, each group of order n/k in groups/order_<m>/ and each
base in bases/, run one `lift` process (exhaust when the orbit-count estimate
is small enough, otherwise sample or skip) and append a JSON line to the
results file.  One process at a time, every process under ulimit -v and timeout.

usage: sweep.py --D 6 --n 134 136 ... [--kmin 1] [--kmax 8] [--exhaust-max 3e8]
                [--sample-seconds 0] [--total-seconds 3000] [--run-timeout 1800]
                [--results results.jsonl] [--bases bases] [--groups groups]
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
VERIFY = os.path.join(HERE, "verify.py")
PY = sys.executable


def read_base(path):
    with open(path) as f:
        k, nslots, ntree = map(int, f.readline().split())
        slots = [tuple(map(int, f.readline().split())) for _ in range(nslots)]
    return k, nslots, ntree, slots


def base_props(path):
    k, nslots, ntree, slots = read_base(path)
    loops = sum(1 for s in slots if s[2] == 1)
    semi = sum(1 for s in slots if s[2] == 2)
    free_group = nslots - ntree - semi
    # bipartite test of the underlying multigraph (type-0 slots only)
    adj = {i: [] for i in range(k)}
    for u, v, t in slots:
        if t == 0:
            adj[u].append(v)
            adj[v].append(u)
    col = {0: 0}
    stack = [0]
    bip = True
    while stack and bip:
        u = stack.pop()
        for w in adj[u]:
            if w not in col:
                col[w] = col[u] ^ 1
                stack.append(w)
            elif col[w] == col[u]:
                bip = False
                break
    return dict(k=k, nslots=nslots, ntree=ntree, loops=loops, semi=semi, free_group=free_group,
                bipartite_plain=(bip and loops == 0 and semi == 0))


def read_group_index(gdir):
    out = []
    with open(os.path.join(gdir, "index.txt")) as f:
        for line in f:
            if line.startswith("#"):
                continue
            m = re.match(r"(g\d+) naut_used=(\d+) (.*?): order (\d+), .*involutions=(\d+)", line)
            if not m:
                continue
            out.append(dict(g=m.group(1), naut=int(m.group(2)), name=m.group(3), m=int(m.group(4)),
                            ninv=int(m.group(5))))
    return out


def load_done(path):
    done = set()
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                try:
                    r = json.loads(line)
                except ValueError:
                    continue
                if r.get("mode") == "skip" and r.get("reason", "").startswith("too big"):
                    continue  # may be rerun with a larger threshold or sampling
                if r.get("mode") == "exhaust" and r.get("exit", 0) != 0:
                    continue  # killed / timed out: rerun
                done.add((r["D"], r["n"], r["k"], r["m"], r["g"], r["base"]))
    return done


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--D", type=int, required=True)
    ap.add_argument("--n", type=int, nargs="+", required=True)
    ap.add_argument("--kmin", type=int, default=1)
    ap.add_argument("--kmax", type=int, default=8)
    ap.add_argument("--exhaust-max", type=float, default=3e8, help="max estimated orbit representatives")
    ap.add_argument("--sample-seconds", type=float, default=0.0)
    ap.add_argument("--total-seconds", type=float, default=3000.0)
    ap.add_argument("--run-timeout", type=float, default=1800.0)
    ap.add_argument("--results", default=os.path.join(HERE, "results.jsonl"))
    ap.add_argument("--bases", default=os.path.join(HERE, "bases"))
    ap.add_argument("--groups", default=os.path.join(HERE, "groups"))
    ap.add_argument("--only-base-re", default=None, help="regex on base name to restrict")
    ap.add_argument("--only-group-re", default=None, help="regex on group name (gN) to restrict")
    ap.add_argument("--force", action="store_true", help="rerun pairs already present in the results file")
    args = ap.parse_args()

    t_start = time.time()
    done = set() if args.force else load_done(args.results)
    soldir = os.path.join(HERE, "solutions")
    os.makedirs(soldir, exist_ok=True)
    lift = os.path.join(HERE, "lift")

    base_names = sorted(x[:-4] for x in os.listdir(args.bases) if re.match(r"^[a-z]\d+_\d+\.txt$", x))
    bases_by_k = {}
    for b in base_names:
        k = int(b[1:].split("_")[0])
        bases_by_k.setdefault(k, []).append(b)

    stats = dict(run=0, skipped=0, exhausted=0, sampled=0, solutions=0, evals=0)
    out = open(args.results, "a")
    stop = False
    for n in args.n:
        for k in range(args.kmin, args.kmax + 1):
            if n % k or k not in bases_by_k:
                continue
            m = n // k
            gdir = os.path.join(args.groups, f"order_{m}")
            if not os.path.isdir(gdir):
                print(f"n={n} k={k}: no groups of order {m} generated", flush=True)
                continue
            groups = read_group_index(gdir)
            for grp in groups:
                if args.only_group_re and not re.search(args.only_group_re, grp["g"]):
                    continue
                gpath = os.path.join(gdir, grp["g"] + ".txt")
                for b in bases_by_k[k]:
                    if args.only_base_re and not re.search(args.only_base_re, b):
                        continue
                    key = (args.D, n, k, m, grp["g"], b)
                    if key in done:
                        continue
                    bp = base_props(os.path.join(args.bases, b + ".txt"))
                    rec = dict(D=args.D, n=n, k=k, m=m, g=grp["g"], gname=grp["name"], naut=grp["naut"],
                               ninv=grp["ninv"], base=b, free_group=bp["free_group"], semi=bp["semi"],
                               loops=bp["loops"])
                    raw = float(m) ** bp["free_group"] * float(grp["ninv"]) ** bp["semi"]
                    est = raw / max(1, grp["naut"])
                    rec["raw"] = raw
                    rec["est"] = est
                    if args.D == 6 and bp["bipartite_plain"]:
                        rec["mode"] = "skip"
                        rec["reason"] = "bipartite lift (<=126 vertices for D=6)"
                    elif bp["semi"] and grp["ninv"] == 0:
                        rec["mode"] = "skip"
                        rec["reason"] = "semi-edges need involutions; group has none"
                        rec["evaluated"] = 0
                    elif est <= args.exhaust_max:
                        rec["mode"] = "exhaust"
                    elif args.sample_seconds > 0:
                        rec["mode"] = "sample"
                    else:
                        rec["mode"] = "skip"
                        rec["reason"] = "too big for exhaust, sampling disabled"
                    if rec["mode"] == "skip":
                        stats["skipped"] += 1
                        out.write(json.dumps(rec) + "\n")
                        out.flush()
                        continue
                    if time.time() - t_start > args.total_seconds:
                        stop = True
                        break
                    prefix = os.path.join(soldir, f"D{args.D}_n{n}_k{k}_m{m}_{grp['g']}_{b}")
                    if rec["mode"] == "exhaust":
                        cmd = [lift, gpath, os.path.join(args.bases, b + ".txt"), str(args.D), "exhaust", "3", prefix]
                    else:
                        cmd = [lift, gpath, os.path.join(args.bases, b + ".txt"), str(args.D), "sample",
                               str(args.sample_seconds), str(stats["run"] + 1), "3", prefix]
                    t0 = time.time()
                    shell = f"ulimit -v 8388608; exec timeout {int(args.run_timeout)} " + " ".join(cmd)
                    p = subprocess.run(["bash", "-c", shell], capture_output=True, text=True)
                    dt = time.time() - t0
                    rec["exit"] = p.returncode
                    rec["wall"] = round(dt, 3)
                    mres = re.search(r"RESULT (.*)", p.stdout)
                    if mres:
                        for kv in mres.group(1).split():
                            kk, vv = kv.split("=")
                            if kk in ("evaluated", "solutions", "best_defect", "nfree", "naut", "ninv"):
                                rec[kk] = int(vv)
                            elif kk == "seconds":
                                rec[kk] = float(vv)
                    else:
                        rec["stderr"] = p.stderr[-300:]
                        rec["stdout_tail"] = p.stdout[-300:]
                    mbest = re.search(r"BEST defect=(\d+) voltages:(.*)", p.stdout)
                    if mbest:
                        rec["best_voltages"] = mbest.group(2).strip()
                    near = re.findall(r"NEAR defect=(\d+) voltages:(.*)", p.stdout)
                    if near:
                        rec["near"] = [[int(d), v.strip()] for d, v in near]
                    sols = re.findall(r"SOLUTION (\S+) voltages:(.*)", p.stdout)
                    if sols:
                        rec["solution_files"] = [s[0] for s in sols]
                        rec["voltages"] = [s[1].strip() for s in sols]
                        stats["solutions"] += len(sols)
                        print(f"*** SOLUTION {rec}", flush=True)
                        for sf, _ in sols:
                            v = subprocess.run([PY, VERIFY, sf, str(args.D)], capture_output=True, text=True)
                            print(v.stdout, flush=True)
                            with open(os.path.join(soldir, "VERIFY.txt"), "a") as vf:
                                vf.write(v.stdout + "\n")
                    stats["run"] += 1
                    stats["evals"] += rec.get("evaluated", 0)
                    if rec["mode"] == "exhaust" and p.returncode == 0:
                        stats["exhausted"] += 1
                    elif rec["mode"] == "sample":
                        stats["sampled"] += 1
                    out.write(json.dumps(rec) + "\n")
                    out.flush()
                if stop:
                    break
                print(f"n={n} k={k} m={m} {grp['g']} {grp['name']}: done; stats {stats} elapsed {time.time() - t_start:.0f}s", flush=True)
            if stop:
                break
        if stop:
            print("TIME BUDGET REACHED, stopping", flush=True)
            break
    out.close()
    print(f"FINAL stats {stats} elapsed {time.time() - t_start:.0f}s", flush=True)


if __name__ == "__main__":
    main()
