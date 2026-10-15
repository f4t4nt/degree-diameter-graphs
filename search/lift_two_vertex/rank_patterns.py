"""Pattern filter and sub-block ranking for the shape-103 family over Z_N x|_r Z_M at diameter D.

usage: rank_patterns.py N M r D            -> table of sub-blocks (i <= j, unit layers): fraction of link-level
                                              patterns (q2, q3) passing R1 (quotient girth >= 4 with distinct link
                                              levels) and R1 & score ratio < 1.01, predicted relative yield
       rank_patterns.py N M r D i j        -> the passing (q2, q3) patterns of sub-block (i, j), one per line
The features are those of pattern_features.py (bi-circulant quotient Q on 2M vertices): quotient girth = shortest
closed non-backtracking walk in Q; score ratio = sum over the 2M cells of N e^(-W_cell/N) divided by the
uniform value 2 n e^(-Moore/n), words of length <= D from both roots.
"""
import math
import sys

import numpy as np

from pattern_features import features


def main():
    N, M, r, D = (int(v) for v in sys.argv[1:5])
    units = [i for i in range(1, (M + 1) // 2) if math.gcd(i, M) == 1]
    if len(sys.argv) > 6:
        i, j = int(sys.argv[5]), int(sys.argv[6])
        for q2 in range(M):
            for q3 in range(M):
                f = features(N, M, D, i, j, [0, q2, q3])
                if f["qgirth"] >= 4 and f["full1"]:
                    print(f"{q2} {q3} score_ratio={f['score_ratio']:.4f} qgirth={f['qgirth']} spheres={f['sph']}")
        return
    rows = []
    for i in units:
        for j in units:
            if j < i:
                continue
            fs = [features(N, M, D, i, j, [0, q2, q3]) for q2 in range(M) for q3 in range(M)]
            r1 = np.mean([f["qgirth"] >= 4 and f["full1"] for f in fs])
            r1s = np.mean([f["qgirth"] >= 4 and f["full1"] and f["score_ratio"] < 1.01 for f in fs])
            dead = np.mean([f["score_ratio"] > 1.2 for f in fs])
            rows.append((i, j, r1, r1s, dead))
    rows.sort(key=lambda t: -t[3])
    print(f"Z_{N} x|_{r} Z_{M}, D = {D}, n = {2 * N * M}: {len(rows)} sub-blocks (i <= j unit layers); "
          f"sorted by the fraction of link patterns passing R1 & score < 1.01")
    print("i j | frac R1 | frac R1 & score<1.01 | frac dead (score>1.2) | note")
    for i, j, r1, r1s, dead in rows:
        note = "equal layers: every lift counted twice (fibre swap)" if i == j else ""
        print(f"{i} {j} | {r1:.3f} | {r1s:.3f} | {dead:.3f} | {note}")


if __name__ == "__main__":
    main()
