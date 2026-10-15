# Degree 5, Diameter 8

A graph with 18,810 vertices, maximum degree 5 and diameter 8. The Moore bound for degree 5 and diameter 8 is 109,226.

## Construction

A regular voltage lift over $`G = \mathbb{Z}_{209} \rtimes_{4} \mathbb{Z}_{45} = \langle x, y \mid x^{209} = y^{45} = 1,\ y x y^{-1} = x^{4} \rangle`$, of order 9,405. The base graph has vertices 0 and 1, with three parallel links $`0 \to 1`$ (voltages $`1`$, $`x y`$, $`x^{136} y^{21}`$), a loop at 0 (voltage $`y^{17}`$) and a loop at 1 (voltage $`x^{18} y^{20}`$).

The lift has vertex set $`\{0, 1\} \times G`$ and joins $`(u, g)`$ to $`(v, g a)`$ for every dart $`u \to v`$ with voltage $`a`$; a loop with voltage $`a`$ joins $`(v, g)`$ to $`(v, g a)`$ and $`(v, g a^{-1})`$. $`G`$ acts on each fibre $`\{v\} \times G`$ by left multiplication.

Vertex $`9405 v + i + 209 j`$ is the pair $`(v, x^i y^j)`$: the first 9,405 vertices lie over base vertex 0 and the other 9,405 over base vertex 1, and vertex $`9405 v`$ is $`(v, e)`$.

## How it was found

Complete $`S_3`$-reduced sub-block sweeps (fixed loop layers) of the shape-103 two-vertex lift family over the best-screened odd-order metacyclic groups; the graph lives over $`\mathbb{Z}_{209} \rtimes_4 \mathbb{Z}_{45}`$. The search code is in [search/lift_two_vertex](../../../search/lift_two_vertex).

## Files

- `deg5-dia8-n18810.edges.txt`: the 47,025 edges, one `u v` per line in the numbering above, grouped by voltage under `#` headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter; the edge list is also checked to be exactly the construction in `construction.json`.

    python3 tools/ddgraph.py check graphs/deg5/dia8
