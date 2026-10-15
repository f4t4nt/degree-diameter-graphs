# Degree 7, Diameter 6

A graph with 12,654 vertices, maximum degree 7 and diameter 6. The Moore bound for degree 7 and diameter 6 is 65,318.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{703} \rtimes_{85} \mathbb{Z}_{18} = \langle x, y \mid x^{703} = y^{18} = 1,\ y x y^{-1} = x^{85} \rangle`$, of order 12,654 (with $`703 = 19 \cdot 37`$), and connection set $`S = \{y^{9}, x y, x^{24} y^{5}, x^{289} y^{8}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`y^{9}`$ is one of the 37 involutions of $`G`$, so the degree is $`1 + 2 \cdot 3 = 7`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 703 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 703`$, $`0 \le j < 18`$), so vertex 0 is the identity and the 18 runs of 703 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 18`$ (see [the numbering figure](../../cayley-numbering.svg)).

## How it was found

Complete sweep of all one-involution connection sets (the involution $`y^9`$ and three inverse pairs) on $`\mathbb{Z}_{703} \rtimes_{85} \mathbb{Z}_{18}`$ up to the automorphisms fixing the involution (exactly one solution). The search code is in [search/cayley_exhaustive](../../../search/cayley_exhaustive).

## Files

- `deg7-dia6-n12654.edges.txt`: the 44,289 edges, one `u v` per line in the numbering above, grouped by generator under `#` headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter; the edge list is also checked to be exactly the construction in `construction.json`.

    python3 tools/ddgraph.py check graphs/deg7/dia6
