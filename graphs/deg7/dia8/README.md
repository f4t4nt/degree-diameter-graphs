# Degree 7, Diameter 8

A graph with 264,440 vertices, maximum degree 7 and diameter 8. The Moore bound for degree 7 and diameter 8 is 2,351,462.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{6611} \rtimes_{193} \mathbb{Z}_{40} = \langle x, y \mid x^{6611} = y^{40} = 1,\ y x y^{-1} = x^{193} \rangle`$, of order 264,440 (with $`6611 = 11 \cdot 601`$), and connection set $`S = \{y^{20}, x y^{4}, x^{4755} y^{9}, x^{2638} y^{6}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`y^{20}`$ is one of the 601 involutions of $`G`$, so the degree is $`1 + 2 \cdot 3 = 7`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 6611 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 6611`$, $`0 \le j < 40`$), so vertex 0 is the identity and the 40 runs of 6,611 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 40`$.

## How it was found

Complete mode-F sweeps of the unit-$`x`$ one-involution Cayley family per (group, level pattern) with the automorphism quotient; $`\mathrm{Cay}(\mathbb{Z}_{6611} \rtimes_{193} \mathbb{Z}_{40})`$ with the involution $`y^{20}`$ and three inverse pairs. The search code is in [search/cayley_exhaustive](../../../search/cayley_exhaustive).

## Files

- `deg7-dia8-n264440.edges.txt.gz`: the 925,540 edges (gzipped), one `u v` per line in the numbering above, grouped by generator under `#` headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The edge list is checked to be exactly the construction in `construction.json`. The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg7/dia8
