# Degree 8, Diameter 9

A graph with 4,471,111 vertices, maximum degree 8 and diameter 9; the Moore bound is 53,804,809.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{66733} \rtimes_{1441} \mathbb{Z}_{67} = \langle x, y \mid x^{66733} = y^{67} = 1,\ y x y^{-1} = x^{1441} \rangle`$, of order 4,471,111 (with 66733 prime), and connection set $`S = \{x^{26424} y^{10}, x^{11037} y^{15}, x^{10559} y^{29}, x^{27701} y^{2}\}`$: each $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for each $`s \in S`$.

$`G`$ has odd order, so it has no involutions, and the degree is $`2 \cdot 4 = 8`$.

Vertex $`i + 66733 j`$ is $`x^i y^j`$ for $`0 \le i < 66733`$ and $`0 \le j < 67`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Groups chosen by the small-complement rule ($`M \approx 0.4\, n^{1/3}`$), simulated annealing with targeted moves on the odd-order group $`\mathbb{Z}_{66733} \rtimes_{1441} \mathbb{Z}_{67}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (17,884,444 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg8/dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg8/dia9
