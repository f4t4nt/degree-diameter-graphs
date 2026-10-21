# Degree 7, Diameter 6

A graph with 12,654 vertices, maximum degree 7 and diameter 6; the Moore bound is 65,318.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{703} \rtimes_{85} \mathbb{Z}_{18} = \langle x, y \mid x^{703} = y^{18} = 1,\ y x y^{-1} = x^{85} \rangle`$, of order 12,654 (with $`703 = 19 \cdot 37`$), and connection set $`S = \{y^{9}, x y, x^{24} y^{5}, x^{289} y^{8}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`y^{9}`$ is one of the 37 involutions of $`G`$, so the degree is $`1 + 2 \cdot 3 = 7`$.

Vertex $`i + 703 j`$ is $`x^i y^j`$ for $`0 \le i < 703`$ and $`0 \le j < 18`$, as in [the numbering figure](../../cayley-numbering.svg).

## Structure

The automorphism group is $`G`$ acting by left multiplication: the graph is a graphical regular representation of $`G`$.

## How it was found

Complete sweep of all one-involution connection sets (the involution $`y^9`$ and three inverse pairs) on $`\mathbb{Z}_{703} \rtimes_{85} \mathbb{Z}_{18}`$ up to the automorphisms fixing the involution (exactly one solution). The search code is in [search/cayley_exhaustive](../../../search/cayley_exhaustive).

## Files

- `deg7-dia6-n12654.edges.txt`: the 44,289 edges as `u v` lines, under a `#` heading per generator.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter; the edge list is also checked against `construction.json`.

    python3 tools/ddgraph.py check graphs/deg7/dia6
