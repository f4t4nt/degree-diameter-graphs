# Degree 5, Diameter 9

A graph with 60,686 vertices, maximum degree 5 and diameter 9. The Moore bound for degree 5 and diameter 9 is 436,906.

## Construction

A regular voltage lift over $`G = \mathbb{Z}_{1597} \rtimes_{3} \mathbb{Z}_{19} = \langle x, y \mid x^{1597} = y^{19} = 1,\ y x y^{-1} = x^{3} \rangle`$, of order 30,343. The base graph has vertices 0 and 1, with three parallel links $`0 \to 1`$ (voltages $`1`$, $`x y^{3}`$, $`x^{1500} y^{9}`$), a loop at 0 (voltage $`y^{2}`$) and a loop at 1 (voltage $`x^{1449} y^{4}`$).

The lift has vertex set $`\{0, 1\} \times G`$ and joins $`(u, g)`$ to $`(v, g a)`$ for every dart $`u \to v`$ with voltage $`a`$; a loop with voltage $`a`$ joins $`(v, g)`$ to $`(v, g a)`$ and $`(v, g a^{-1})`$. $`G`$ acts on each fibre $`\{v\} \times G`$ by left multiplication.

Vertex $`30343 v + i + 1597 j`$ is the pair $`(v, x^i y^j)`$: the first 30,343 vertices lie over base vertex 0 and the other 30,343 over base vertex 1, and vertex $`30343 v`$ is $`(v, e)`$.

## How it was found

Complete sweep of the S_3-reduced sub-block with loop layers 2 and 4 of the shape-103 two-vertex lift family over the odd-order group Z_1597 x|_3 Z_19. The search code is in [search/lift_two_vertex](../../../search/lift_two_vertex).

## Files

- `deg5-dia9-n60686.edges.txt`: the 151,715 edges, one `u v` per line in the numbering above, grouped by voltage under `#` headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter; the edge list is also checked to be exactly the construction in `construction.json`.

    python3 tools/ddgraph.py check graphs/deg5/dia9
