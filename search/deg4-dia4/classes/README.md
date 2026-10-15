# The 36 Graphs Of Degree 4 And Diameter 4 On 104 Vertices

I enumerated every diameter-4 realisation of the three tight three-block skeletons that admit one
(see `../proofs/TEMPLATE_CALCULUS.md`, sections 8 and 9) and sorted the results into isomorphism
classes with nauty: 8 + 18 + 10 = 36 classes, which this directory lists one file per class.

- `class00.edges.txt` to `class35.edges.txt`: one representative per class, 104 vertices, 208 edges,
  one `u v` pair per line, `#` comment lines grouping the edges by type. Classes 0 to 34 are the
  graphs the search campaign had already found, in the order of its corpus; class 35 is the one
  that only the enumeration produced. Classes 29 and 35 have automorphism group of order 32 and
  are the two pure templates (every connector adjacent to one port per block and nothing bent);
  class 20 has order 2; the other 33 are asymmetric. The graph the campaign reported first is
  class 10, its sibling repair is class 11, the independent second seed gave class 27, and the
  graph in `graphs/deg4/dia4/` of this repository is class 29.
- Numbering: block $`A`$ is $`0, \ldots, 23`$, $`B`$ is $`24, \ldots, 47`$, $`C`$ is $`48, \ldots, 71`$.
  Inside a block the hubs are $`L_0, \ldots, L_3`$ at offsets $`0, \ldots, 3`$ and $`R_0, \ldots, R_3`$
  at offsets $`4, \ldots, 7`$, and the port between $`L_x`$ and $`R_y`$ is at offset $`8 + 4x + y`$
  (so the 32 hub-port edges of a block are forced by the numbering). The connectors are
  $`72, \ldots, 103`$ in matched pairs: pair $`k`$ is the adjacent pair $`(72 + 2k, 73 + 2k)`$. The
  remaining content of a class is its 96 port-connector edges.
- `export_classes.py`: produced the files from the enumeration output of
  `../sat_closure/classify104.py` (which labels vertices as hub, port or connector) by the map above.
- `check_classes.py`: plain BFS (order, degrees, diameter, pairs beyond 4) and pynauty certificates
  for a set of edge files, reporting any two that are isomorphic; `--match FILE... --` also says
  which class other edge files belong to.

Verification done on these files: every file has 104 vertices, is 4-regular, simple, connected,
of diameter 4 with 0 pairs beyond distance 4 (both the Python BFS in `check_classes.py` and the C
verifier `../lift_repair/verify_bfs.c`), the 36 certificates are pairwise distinct, and each file
is isomorphic to the campaign's representative of the same index.

To run it yourself:

```
python check_classes.py class*.edges.txt
cc -O2 -o verify_bfs ../lift_repair/verify_bfs.c && for f in class*.edges.txt; do ./verify_bfs $f 4 | tail -1; done
python check_classes.py --match some_other_104.edges -- class*.edges.txt   # which class is it?
```
