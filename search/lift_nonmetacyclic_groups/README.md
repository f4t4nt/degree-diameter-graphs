# Lift Nonmetacyclic Groups

I run the shape-103 two-vertex lift sweep and the one-involution Cayley sweep over groups $`G = K \rtimes C_M`$ whose
odd-order kernel $`K`$ is a direct product of blocks (cyclic, elementary abelian $`\mathbb{F}_p^e`$ with a
$`\mathrm{GL}`$ action, Frobenius $`\mathbb{F}_p^2 \rtimes C_k`$, extraspecial $`p^{1+2}_+`$,
$`\mathbb{Z}_{p^a} \times \mathbb{Z}_{p^b}`$), not only the metacyclic ones, selecting the groups of a window by their
conjugacy-class count $`k(G)`$ (the lower the commuting probability $`k(G)/|G|`$, the softer the family). The complete
shape-103 sweep of $`((\mathbb{F}_{25} \rtimes C_3) \times \mathbb{Z}_5) \rtimes C_8`$, the order-3000 group with the
lowest $`k(G)`$ of its window, gave the (5,7) graph.

Files:

- `gf.py`: builds a group from a JSON block spec (multiplication table of $`K`$, the automorphism $`\alpha`$ of $`y`$)
  into a `.grp` file and reports $`|Z(G)|`$, involutions, element orders; `gf.py info FILE.grp` reads one back.
- `enum_family.py`, `enum_groups.py`: every block decomposition and automorphism class of a window
  $`(\mathrm{LO}, \mathrm{HI}]`$ with the filters of the two families (`--cell cay7`: $`M`$ even, one involution class;
  `--cell lift`: any $`M`$), as JSONL specs.
- `aut_k.py`: automorphisms of $`K`$ commuting with $`\alpha`$ (units of cyclic blocks, $`\mathrm{GL}`$ centralisers) as
  an `.aut` file that shrinks the sweep; any subgroup of $`\mathrm{Aut}(G)`$ gives a complete reduction.
- `table_sweep.c`: the sweeper on the full right-multiplication table (groups up to about 20000 elements): `--family
  lift` (shape 103, reduced by $`\mathrm{Inn}(G) \circ B`$, loop inversions, the $`S_3`$ of link re-gaugings and the
  vertex swap) and `--family cay7` ($`\mathrm{Cay}(G, \{y^{M/2}, a, b, c\}^{\pm})`$); modes `info`, `random`,
  `exhaustive` (with `--chunk i N`), `eval`, `edges` and `cover` (brute-force completeness on small groups).
- `tablefree_sweep.c`, `tfgroup.py`: the same two families without a table, for $`K = S \times \mathbb{Z}_q`$ with the
  cyclic part fused by CRT and a bit-parallel BFS over rows of $`q`$ bits (large groups); `tfgroup.py` writes its `.tfg`
  input.
- `make_edges.py`: independent edge-list rebuild of a FOUND set from the `.grp` file's spec.

Element $`(k, m) = k y^m`$ has index $`k M + m`$ ($`k`$ the mixed-radix block index in spec order); a lift vertex
$`(f, g)`$ has index $`f\,|G| + \mathrm{index}(g)`$. Output lines: `FOUND diam<=D (defect 0) family=.. set=[..]`,
`NEAR`, `PROGRESS`, `EXHAUSTIVE`.

To run it yourself (rebuilds the (5,7) group, re-evaluates the (5,7) lift exactly and brute-force checks the
completeness of the reduction on a small group; all in seconds):

```
pip install -r requirements.txt && python3 gf.py build '{"M":8,"blocks":[{"t":"frob","p":5,"B":[[0,4],[1,4]],"A":[[1,1],[2,4]],"e":2},{"t":"cyc","q":5,"u":2}]}' g3000.grp
# expected last line: {"n": 3000, "nK": 375, "M": 8, "Z": 1, "inv": 25, "ord_alpha": 8}
gcc -O3 -march=native -o table_sweep table_sweep.c -lm && ./table_sweep g3000.grp --family lift --D 7 --mode eval --set 1 316 172 1159
# expected: EVAL family=lift set=[1,316,172,1159] D=7 root0 defect=0 profile=1 5 20 80 319 1215 2885 1475 | root1 defect=0 profile=1 5 20 80 320 1212 2900 1462 | total defect=0
python3 gf.py build '{"M":4,"blocks":[{"t":"vec","p":3,"A":[[0,1],[2,0]]}]}' g36.grp && ./table_sweep g36.grp --family lift --D 4 --mode cover
# expected: COVER total=804440 uncovered=0 OK
```

The complete (5,7) sweep is `python3 aut_k.py '<the spec above>' g3000.aut` followed by `./table_sweep g3000.grp --aut
g3000.aut --family lift --D 7 --mode exhaustive --chunk i 15 --near 4` for $`i = 0, \ldots, 14`$ (about
$`1.9 \times 10^8`$ candidates in total).

## Graphs found

- [(5, 7)](../../graphs/deg5/dia7): 6,000 vertices
