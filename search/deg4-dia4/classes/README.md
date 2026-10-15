# The 36 Graphs Of Degree 4 And Diameter 4 On 104 Vertices

I enumerate every diameter-4 realisation of the three tight three-block skeletons that admit one
(see `../proofs/TEMPLATE_CALCULUS.md`, sections 8 and 9) and sort the results into isomorphism
classes with nauty: 8 + 18 + 10 = 36 classes, one file per class in this directory.

- `class00.edges.txt` to `class35.edges.txt`: one representative per class, 104 vertices, 208 edges,
  one `u v` pair per line, `#` comment lines grouping the edges by type. The lift-and-repair search
  (`../lift_repair/`) finds classes 0 to 34; class 35 appears only in the enumeration. Classes 29
  and 35 have automorphism group of order 32 and are the two pure templates (every connector
  adjacent to one port per block and nothing bent); class 20 has order 2; the other 33 are
  asymmetric. The graph in [`graphs/deg4/dia4/`](../../../graphs/deg4/dia4) of this repository is
  class 29.
- Numbering: block $`A`$ is $`0, \ldots, 23`$, $`B`$ is $`24, \ldots, 47`$, $`C`$ is $`48, \ldots, 71`$.
  Inside a block the hubs are $`L_0, \ldots, L_3`$ at offsets $`0, \ldots, 3`$ and $`R_0, \ldots, R_3`$
  at offsets $`4, \ldots, 7`$, and the port between $`L_x`$ and $`R_y`$ is at offset $`8 + 4x + y`$
  (so the numbering forces the 32 hub-port edges of a block). The connectors are
  $`72, \ldots, 103`$ in matched pairs: pair $`k`$ is the adjacent pair $`(72 + 2k, 73 + 2k)`$. The
  rest of a class is its 96 port-connector edges.
- `export_classes.py`: writes the files from the enumeration output of
  `../sat_closure/classify104.py` (which labels vertices as hub, port or connector) by the map above.
- `check_classes.py`: plain BFS (order, degrees, diameter, pairs beyond 4) and pynauty certificates
  for a set of edge files, reporting isomorphic pairs; `--match FILE... --` also names the class of
  other edge files.

Verification of these files: every file has 104 vertices, is 4-regular, simple, connected, of
diameter 4 with 0 pairs beyond distance 4 (both the Python BFS in `check_classes.py` and the
repository's C verifier `tools/verify/verify_bfs.c`), and the 36 certificates are pairwise
distinct.

To run it yourself:

```
python check_classes.py class*.edges.txt
cc -O2 -o verify_bfs ../../../tools/verify/verify_bfs.c && for f in class*.edges.txt; do ./verify_bfs $f 4 | tail -1; done
python check_classes.py --match some_other_104.edges -- class*.edges.txt   # which class is it?
```
