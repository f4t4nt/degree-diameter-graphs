# Degree 8, Diameter 10

A graph with 26,864,290 vertices, maximum degree 8 and diameter 10. The Moore bound for degree 8 and diameter 10 is 376,633,665.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{141391} \rtimes_{240} \mathbb{Z}_{190} = \langle x, y \mid x^{141391} = y^{190} = 1,\ y x y^{-1} = x^{240} \rangle`$, of order 26,864,290 (with $`141391 = 31 \cdot 4561`$), and connection set $`S = \{x^{50270} y^{33}, x^{38396} y^{111}, x^{55026} y^{57}, x^{60295} y^{14}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`S`$ contains no involution ($`G`$ has 31 involutions), so the degree is $`2 \cdot 4 = 8`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 141391 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 141391`$, $`0 \le j < 190`$), so vertex 0 is the identity and the 190 runs of 141,391 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 190`$ (see [the numbering figure](../../cayley-numbering.svg)).

## How it was found

Even-$`M`$ group with a tiny involution class, coincidence-free level pattern, simulated annealing with targeted moves on $`\mathbb{Z}_{141391} \rtimes_{240} \mathbb{Z}_{190}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (107,457,160 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg8/dia10` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg8/dia10
