"""Random screen of whole shape-103 families over several groups Z_N x|_r Z_M with one statistic: the u_D(root 0)
histogram of SAMPLES random candidates whose loops a, b lie in "unit layers" (y-exponent coprime to M; the other
layers are uniformly worse) and whose links c2, c3 are uniform over G.  Compiles the bitlift binary
(../lift_two_vertex/bitlift.c) for the right NW = ceil(N/64) on demand.

usage: screen_groups.py D outprefix [LINES RUNB] < groups  (lines "n N M r ...", '#' comments ignored)
"""
import math
import os
import random
import re
import subprocess
import sys

BASE = "k=2;deg=5;loop:0;link*:0,1;link:0,1;link:0,1;loop:1"


def binary_for(N):
    nw = (N + 63) // 64
    fn = f"./bitlift_nw{nw}"
    if not os.path.exists(fn):
        src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lift_two_vertex", "bitlift.c")
        subprocess.run(["gcc", "-O3", "-march=native", f"-DNW={nw}", "-o", fn, src, "-lm"], check=True)
    return fn


def main():
    D = int(sys.argv[1])
    outprefix = sys.argv[2]
    lines = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
    runb = int(sys.argv[4]) if len(sys.argv) > 4 else 50
    print(f"# group screen D={D} lines={lines} runb={runb} per group ({lines * runb} candidates), loops in unit layers")
    print("# n N M r : simple best rate  hist(u_D root0: [0] [1] [2-3] [4-7] [8-15] [16-31] [32-63]) rest(>=64)")
    sys.stdout.flush()
    for ln in sys.stdin:
        if ln.startswith("#") or not ln.strip():
            continue
        n, N, M, r = (int(v) for v in ln.split()[:4])
        m = N * M
        rng = random.Random(1000003 * N + 17 * M + r)
        ulayers = [i for i in range(1, M) if math.gcd(i, M) == 1]
        fn = f"{outprefix}_{N}_{M}_{r}.txt"
        with open(fn, "w") as f:
            for _ in range(lines):
                a = rng.randrange(N) + N * rng.choice(ulayers)
                c2 = rng.randrange(1, m)
                c3 = rng.randrange(1, m)
                while c3 == c2:
                    c3 = rng.randrange(1, m)
                lo = rng.randrange(N - runb) + N * rng.choice(ulayers)
                f.write(f"{a} {c2} {c3} {lo} {lo + runb}\n")
        binary = binary_for(N)
        cmd = [binary, "--meta", str(N), str(M), str(r), "--base", BASE, "--D", str(D), "--mode", "prefix", "--prefix", fn, "--report", "2"]
        res = subprocess.run(cmd, capture_output=True, text=True)
        out = res.stdout
        mt = re.search(r"candidates=(\d+) simple=(\d+) .*rate=(\d+) cand/s best_defect=(\d+)", out)
        if not mt:
            print(f"{n} {N} {M} {r} : FAILED {res.stderr.strip()[:200]}")
            sys.stdout.flush()
            continue
        hist = {}
        for b, c in re.findall(r"\[(?:<2\^)?(\d+)\]=(\d+)", out):
            hist[int(b)] = int(c)
        bins = [hist.get(b, 0) for b in range(7)]
        rest = sum(v for b, v in hist.items() if b >= 7)
        near = [x for x in out.splitlines() if x.startswith("NEAR") or x.startswith("FOUND")]
        print(f"{n} {N} {M} {r} : {mt.group(2)} best={mt.group(4)} rate={mt.group(3)}  {bins} rest={rest}" + (f"  {near}" if near else ""))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
