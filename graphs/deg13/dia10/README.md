# Degree 13, Diameter 10

A graph with 3,789,298,110 vertices, maximum degree 13 and diameter 10; the Moore bound is 73,175,066,810.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{9716149} \rtimes_{514348} \mathbb{Z}_{390} = \langle x, y \mid x^{9716149} = y^{390} = 1,\ y x y^{-1} = x^{514348} \rangle`$, of order 3,789,298,110 (with $`9716149 = 859 \cdot 11311`$), and connection set $`S = \{x^{3145658} y^{195}, x^{4424168} y^{343}, x^{5603880} y^{4}, x^{9317141} y^{59}, x^{2851313} y, x^{373050} y^{31}, x^{4695647} y^{146}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`x^{3145658} y^{195}`$ is one of the 11,311 involutions of $`G`$, so the degree is $`1 + 2 \cdot 6 = 13`$.

Vertex $`i + 9716149 j`$ is $`x^i y^j`$ for $`0 \le i < 9716149`$ and $`0 \le j < 390`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Two-plane bitmap simulated annealing with targeted moves and a movable involution in the one-involution Cayley family on $`\mathbb{Z}_{9716149} \rtimes_{514348} \mathbb{Z}_{390}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (24,630,437,715 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg13/dia10` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search, which `check` skips unless asked:

    python3 tools/ddgraph.py check graphs/deg13/dia10 --max-vertices 3789298110
