# Degree 8, Diameter 9

A graph with 4,471,111 vertices, maximum degree 8 and diameter 9. The Moore bound for degree 8 and diameter 9 is 53,804,809.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{66733} \rtimes_{1441} \mathbb{Z}_{67} = \langle x, y \mid x^{66733} = y^{67} = 1,\ y x y^{-1} = x^{1441} \rangle`$, of order 4,471,111 (with 66733 prime), and connection set $`S = \{x^{26424} y^{10}, x^{11037} y^{15}, x^{10559} y^{29}, x^{27701} y^{2}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`G`$ has odd order, so it has no involutions, and the degree is $`2 \cdot 4 = 8`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 66733 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 66733`$, $`0 \le j < 67`$), so vertex 0 is the identity and the 67 runs of 66,733 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 67`$.

## How it was found

Groups chosen by the small-complement rule (M about 0.4 n^(1/3)), simulated annealing with targeted moves on the odd-order group Z_66733 x|_1441 Z_67. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (17,884,444 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg8/dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg8/dia9
