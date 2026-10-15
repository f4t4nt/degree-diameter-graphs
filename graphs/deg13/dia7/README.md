# Degree 13, Diameter 7

A graph with 3,691,038 vertices, maximum degree 13 and diameter 7; the Moore bound is 42,346,682.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{47321} \rtimes_{3981} \mathbb{Z}_{78} = \langle x, y \mid x^{47321} = y^{78} = 1,\ y x y^{-1} = x^{3981} \rangle`$, of order 3,691,038 (with $`47321 = 79 \cdot 599`$), and connection set $`S = \{x^{46136} y^{39}, x^{28946} y^{31}, x^{23110} y^{28}, x^{30322} y^{32}, x^{27810} y^{11}, x^{25477} y^{37}, x^{46989} y^{12}\}`$: each $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for each $`s \in S`$.

$`x^{46136} y^{39}`$ is one of the 599 involutions of $`G`$, so the degree is $`1 + 2 \cdot 6 = 13`$.

Vertex $`i + 47321 j`$ is $`x^i y^j`$ for $`0 \le i < 47321`$ and $`0 \le j < 78`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Simulated annealing with targeted moves and a movable involution in the one-involution Cayley family on $`\mathbb{Z}_{47321} \rtimes_{3981} \mathbb{Z}_{78}`$, the group chosen by the lower-quartile softness screen. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (23,991,747 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg13/dia7` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg13/dia7
