# Degree 11, Diameter 7

A graph with 1,197,748 vertices, maximum degree 11 and diameter 7; the Moore bound is 12,222,222.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{13019} \rtimes_{770} \mathbb{Z}_{92} = \langle x, y \mid x^{13019} = y^{92} = 1,\ y x y^{-1} = x^{770} \rangle`$, of order 1,197,748 (with $`13019 = 47 \cdot 277`$), and connection set $`S = \{x^{7614} y^{46}, x^{4646} y^{8}, x^{12909} y^{15}, x^{7519} y^{6}, x^{4020} y^{18}, x^{5180} y^{57}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`x^{7614} y^{46}`$ is one of the 277 involutions of $`G`$, so the degree is $`1 + 2 \cdot 5 = 11`$.

Vertex $`i + 13019 j`$ is $`x^i y^j`$ for $`0 \le i < 13019`$ and $`0 \le j < 92`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Simulated annealing with level patterns, targeted moves and a movable involution in the one-involution Cayley family on $`\mathbb{Z}_{13019} \rtimes_{770} \mathbb{Z}_{92}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (6,587,614 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg11/dia7` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg11/dia7
