# Cayley Exhaustive

I sweep complete families of Cayley graphs on metacyclic groups Z_N x|_r Z_M, one representative per orbit of the
group's automorphisms, so that a negative result closes the family and a hit is found wherever it sits: the
one-involution odd-degree family (every set {t, K inverse pairs} up to the automorphisms fixing t), the even-degree
3-pair family (level patterns, every x-assignment) and the fully reduced 3-pair list of a single group. These sweeps
produced the (6,8), (6,10), (7,6), (7,8) and (9,5) records.

Files:

- `sweep_odd.c`: complete sweep of Cay(G, {y^(M/2)} u {a_1^+-1, ..., a_K^+-1}) for N odd and M even, with the
  stabiliser-chain lexicographic-minimum test, chunking (`--c0 lo hi`, `--c1 lo hi`), `--count`, a per-class softness
  screen and brute-force `--check` on small groups.
- `lift.c`: bit-parallel evaluator with an exhaustive prefix mode; a `Y la lb lc` line sweeps every x-assignment of a
  level pattern of three inverse pairs (N + 2 candidates for la > 0, N for la = 0); also random, sa and completion modes.
- `make_ylist.py`: every admissible level pattern of a group, in chunk files with a manifest of exact candidate counts
  (the completeness certificate: lift's `EXHAUSTIVE ... candidates=` totals must equal it).
- `gen_cay.py`, `check_cay.py`: fully symmetry-reduced `E a b L` / `T a b c` candidate lists for 3-pair Cayley graphs on
  one group (automorphisms, inversions, generator permutations) for `lift --mode exhaustive --prefix`; brute-force check.
- `gen_edges.py`: streams the edge list of a Cayley graph or 2-vertex lift given N M r and the generators. Mode F of
  `../cayley_metacyclic_anneal/cayley_sa.c` (`--s2all --s2bound`, quotient `--quot`) is the complete sweep of a
  one-involution family with a fixed level pattern; it is built there.

To run it yourself (the controls re-find the (7,6) graph on 12654 vertices, the unique solution of its group, in about
20 seconds; the known (6,9) graph on 331387 vertices from its level pattern; the known (7,7) graph on 56952 vertices
with the mode-F quotient in a few seconds):

```
gcc -O3 -o sweep_odd sweep_odd.c -lm && ./sweep_odd 703 18 85 3 6
# expected: FOUND 0 : 1 1 24 5 289 8   then   RESULT N=703 M=18 s=85 ... n=12654 K=3 D=6 ... candidates=4364198 min_defect=0 n_min=1 found0=1
gcc -O3 -march=native -o lift lift.c -lm && printf 'Y 2 18 24\n' > ypat.txt
./lift --meta 6763 49 41 --base "k=1;loop:0;loop:0;loop:0" --deg 6 --D 9 --mode exhaustive --prefix ypat.txt --near 4
# expected: FOUND diam<=9 (defect 0) group=Z6763 x|_41 Z49 ... voltages=[13526,121735,165316]   then   EXHAUSTIVE ... candidates=6765 ... found=1 complete=1
../cayley_metacyclic_anneal/cayley_sa 2373 24 44 7 3 0 --invol --xunitI --tx 0 --init 1 1 2 2 4 7 --sweep2 800 F --J2 2 --pairs 1-2 --s2bound 8 --prerow_max 8 --s2all --seed 7 --quot
# expected: FOUND N=2373 M=24 r=44 D=7 k=3 c=0 degree 7 n=56952 + involution (last pair) : 1 1 120 2 728 7 0 12
#           SWEEP2 F done: 385666 candidates, 0 ties, 1 improved ...
```

Voltages in `lift` lines are element indices i + N j of x^i y^j; `make_ylist.py 6763 49 - ylist 2000000` writes all 2884
patterns of that group in chunks. Compile `sweep_odd` without `-march=native` (a known slow-code interaction).
