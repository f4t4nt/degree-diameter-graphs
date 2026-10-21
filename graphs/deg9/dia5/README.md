# Degree 9, Diameter 5

A graph with 9,020 vertices, maximum degree 9 and diameter 5; the Moore bound is 42,130.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{451} \rtimes_{46} \mathbb{Z}_{20} = \langle x, y \mid x^{451} = y^{20} = 1,\ y x y^{-1} = x^{46} \rangle`$, of order 9,020 (with $`451 = 11 \cdot 41`$), and the connection set $`S = \{y^{10}, x y, x^{316} y^{6}, x^{13} y^{7}, x^{426} y^{9}\}`$: every $`g \in G`$ is adjacent to $`g s`$ and $`g s^{-1}`$ for every $`s \in S`$.

$`y^{10}`$ is one of the 41 involutions of $`G`$, so the degree is $`1 + 2 \cdot 4 = 9`$.

Vertex $`i + 451 j`$ is $`x^i y^j`$ for $`0 \le i < 451`$ and $`0 \le j < 20`$, as in [the numbering figure](../../cayley-numbering.svg).

## Structure

The automorphism group is $`G`$ acting by left multiplication: the graph is a graphical regular representation of $`G`$.

## How it was found

Complete chunked sweep of the one-involution Cayley family (the involution $`y^{10}`$ and four inverse pairs) on the groups $`\mathbb{Z}_{451} \rtimes_s \mathbb{Z}_{20}`$ up to automorphisms (two inequivalent solutions on two twists $`s`$). The search code is in [search/cayley_exhaustive](../../../search/cayley_exhaustive).

## Files

- `deg9-dia5-n9020.edges.txt`: the 40,590 edges as `u v` lines, under a `#` heading per generator.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter; the edge list is also checked against `construction.json`.

    python3 tools/ddgraph.py check graphs/deg9/dia5
