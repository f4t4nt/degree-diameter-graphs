# Degree 3, Diameter 7

A graph with 200 vertices, maximum degree 3 and diameter 7. The Moore bound for degree 3 and diameter 7 is 382.

## Construction

A regular voltage lift over $`G = \mathbb{Z}_{5} \rtimes_{4} \mathbb{Z}_{8} = \langle x, y \mid x^{5} = y^{8} = 1,\ y x y^{-1} = x^{4} \rangle`$, of order 40. The base graph has vertices 0, 1, 2, 3 and 4, with links $`0 \to 1`$, $`0 \to 2`$, $`0 \to 3`$ and $`1 \to 4`$ (voltage $`1`$), a link $`1 \to 2`$ (voltage $`y`$), a semi-edge at 2 (voltage $`y^{4}`$, an involution), a loop at 3 (voltage $`x y^{2}`$) and a loop at 4 (voltage $`x^{2} y^{2}`$).

The lift has vertex set $`\{0, 1, 2, 3, 4\} \times G`$ and joins $`(u, g)`$ to $`(v, g a)`$ for every dart $`u \to v`$ with voltage $`a`$; a loop with voltage $`a`$ joins $`(v, g)`$ to $`(v, g a)`$ and $`(v, g a^{-1})`$, and a semi-edge with involution $`t`$ joins $`(v, g)`$ to $`(v, g t)`$. $`G`$ acts on each fibre $`\{v\} \times G`$ by left multiplication.

Vertex $`40 v + i + 5 j`$ is the pair $`(v, x^i y^j)`$: each run of 40 consecutive numbers lies over one base vertex, in order, and vertex $`40 v`$ is $`(v, e)`$.

## How it was found

Exhaustive sweep of the voltage assignments of every 5-vertex cubic base graph (links, loops and semi-edges) over the groups of order 40; the graph is the lift of the base with edges 0-1, 0-2, 0-3, 1-4, 1-2, a semi-edge at 2 and loops at 3 and 4 over $`\mathbb{Z}_5 \rtimes \mathbb{Z}_8`$. The search code is in [search/cubic_lifts](../../../search/cubic_lifts).

## Files

- `deg3-dia7-n200.edges.txt`: the 300 edges, one `u v` per line in the numbering above, grouped by voltage under `#` headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter; the edge list is also checked to be exactly the construction in `construction.json`.

    python3 tools/ddgraph.py check graphs/deg3/dia7
