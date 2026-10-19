# Degree 10, Diameter 9

A graph with 34,913,494 vertices, maximum degree 10 and diameter 9; the Moore bound is 484,275,611.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{226711} \rtimes_{8084} \mathbb{Z}_{154} = \langle x, y \mid x^{226711} = y^{154} = 1,\ y x y^{-1} = x^{8084} \rangle`$, of order 34,913,494 (with $`226711 = 23 \cdot 9857`$), and connection set $`S = \{x^{17400} y^{45}, x^{48393} y^{74}, x^{206671} y, x^{42957} y^{145}, x^{57429} y^{87}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`S`$ contains no involution ($`G`$ has 23 involutions), so the degree is $`2 \cdot 5 = 10`$.

Vertex $`i + 226711 j`$ is $`x^i y^j`$ for $`0 \le i < 226711`$ and $`0 \le j < 154`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Even-$`M`$ group with a tiny involution class chosen by a random-set softness screen of a window of group orders, coincidence-free level pattern, simulated annealing with targeted moves on $`\mathbb{Z}_{226711} \rtimes_{8084} \mathbb{Z}_{154}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (174,567,470 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg10/dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg10/dia9
