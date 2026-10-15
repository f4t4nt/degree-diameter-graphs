# Degree 9, Diameter 9

A graph with 13,219,318 vertices, maximum degree 9 and diameter 9; the Moore bound is 172,565,650.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{134891} \rtimes_{4746} \mathbb{Z}_{98} = \langle x, y \mid x^{134891} = y^{98} = 1,\ y x y^{-1} = x^{4746} \rangle`$, of order 13,219,318 (with $`134891 = 43 \cdot 3137`$), and connection set $`S = \{x^{48203} y^{49}, x^{74247} y, x^{51808} y^{79}, x^{121079} y^{27}, x^{71015} y^{82}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`x^{48203} y^{49}`$ is one of the 3,137 involutions of $`G`$, so the degree is $`1 + 2 \cdot 4 = 9`$.

Vertex $`i + 134891 j`$ is $`x^i y^j`$ for $`0 \le i < 134891`$ and $`0 \le j < 98`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Simulated annealing with targeted moves and a movable involution in the one-involution Cayley family on $`\mathbb{Z}_{134891} \rtimes_{4746} \mathbb{Z}_{98}`$, the group chosen by the lower-quartile softness screen. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (59,486,931 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg9/dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg9/dia9
