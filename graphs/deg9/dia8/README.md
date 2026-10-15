# Degree 9, Diameter 8

A graph with 1,961,370 vertices, maximum degree 9 and diameter 8; the Moore bound is 21,570,706.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{21793} \rtimes_{34} \mathbb{Z}_{90} = \langle x, y \mid x^{21793} = y^{90} = 1,\ y x y^{-1} = x^{34} \rangle`$, of order 1,961,370 (with $`21793 = 19 \cdot 31 \cdot 37`$), and connection set $`S = \{x^{18500} y^{45}, x^{10329} y^{8}, x^{6476} y^{9}, x^{17755} y^{5}, x^{3186} y^{39}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`x^{18500} y^{45}`$ is one of the 589 involutions of $`G`$, so the degree is $`1 + 2 \cdot 4 = 9`$.

Vertex $`i + 21793 j`$ is $`x^i y^j`$ for $`0 \le i < 21793`$ and $`0 \le j < 90`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Simulated annealing with level patterns, targeted moves and a movable involution in the one-involution Cayley family on $`\mathbb{Z}_{21793} \rtimes_{34} \mathbb{Z}_{90}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (8,826,165 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg9/dia8` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg9/dia8
