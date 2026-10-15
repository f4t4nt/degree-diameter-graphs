# Degree 4, Diameter 6

A graph with 760 vertices, maximum degree 4 and diameter 6. The Moore bound for degree 4 and diameter 6 is 1,457.

## Construction

The incidence graph of the generalized hexagon $`H(3)`$ (364 points and 364 lines; 4-regular, bipartite, girth 12, diameter 6) with 32 flags subdivided: a new vertex is placed on each of 32 point-line incidence edges and joined to that point and that line, and the 32 new vertices carry a 2-regular graph $`F`$ with cycles of lengths 8, 8 and 16.

The 32 flags are the union of two 16-flag orbits of a Sylow 2-subgroup (order 128) of $`\mathrm{Aut}(H(3)) = G_2(3){:}2`$ inside the centraliser of an involution. $`F`$ has no edge between the two orbits, so the graph is two 744-vertex solutions glued together.

Vertices 0 to 363 are the points of $`H(3)`$, 364 to 727 its lines and 728 to 759 the 32 new vertices.

## How it was found

Incidence graph of the generalized hexagon $`H(3)`$ (728 vertices) with 32 flags subdivided and a 2-regular graph on the new vertices; the flag set is the union of two 16-flag orbits of a Sylow 2-subgroup of the centraliser of an involution of $`G_2(3)`$, and the 2-factor is found by exact-cost simulated annealing with the flag set fixed. The search code is in [search/hexagon_insertion](../../../search/hexagon_insertion).

## Files

- `deg4-dia6-n760.edges.txt`: the 1,520 edges, one `u v` per line in the numbering above, grouped by edge type under `#` headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter.

    python3 tools/ddgraph.py check graphs/deg4/dia6
