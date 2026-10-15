# Degree 7, Diameter 8

A graph with 264,440 vertices, maximum degree 7 and diameter 8; the Moore bound is 2,351,462.

## Construction

The Cayley graph of the group $`G = \mathbb{Z}_{6611} \rtimes_{193} \mathbb{Z}_{40} = \langle x, y \mid x^{6611} = y^{40} = 1,\ y x y^{-1} = x^{193} \rangle`$, of order 264,440 (with $`6611 = 11 \cdot 601`$), and the connection set $`S = \{y^{20}, x y^{4}, x^{4755} y^{9}, x^{2638} y^{6}\}`$: each $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for each $`s \in S`$.

$`y^{20}`$ is one of the 601 involutions of $`G`$, so the degree is $`1 + 2 \cdot 3 = 7`$.

Vertex $`i + 6611 j`$ is $`x^i y^j`$ for $`0 \le i < 6611`$ and $`0 \le j < 40`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Complete mode-F sweeps of the unit-$`x`$ one-involution Cayley family per (group, level pattern) with the automorphism quotient; $`\mathrm{Cay}(\mathbb{Z}_{6611} \rtimes_{193} \mathbb{Z}_{40})`$ with the involution $`y^{20}`$ and three inverse pairs. The search code is in [search/cayley_exhaustive](../../../search/cayley_exhaustive).

## Files

- `deg7-dia8-n264440.edges.txt.gz`: the 925,540 edges as `u v` lines, under a `#` heading per generator.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The edge list is checked against `construction.json`. The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg7/dia8
