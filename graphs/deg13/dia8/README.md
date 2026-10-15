# Degree 13, Diameter 8

A graph with 36,275,160 vertices, maximum degree 13 and diameter 8; the Moore bound is 508,160,186.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{302293} \rtimes_{7818} \mathbb{Z}_{120} = \langle x, y \mid x^{302293} = y^{120} = 1,\ y x y^{-1} = x^{7818} \rangle`$, of order 36,275,160 (with $`302293 = 41 \cdot 73 \cdot 101`$), and connection set $`S = \{x^{164125} y^{60}, x^{129522} y^{3}, x^{47580} y^{27}, x^{218857} y^{25}, x^{34159} y^{56}, x^{366} y^{18}, x^{71547} y^{44}\}`$: each $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for each $`s \in S`$.

$`x^{164125} y^{60}`$ is one of the 2,993 involutions of $`G`$, so the degree is $`1 + 2 \cdot 6 = 13`$.

Vertex $`i + 302293 j`$ is $`x^i y^j`$ for $`0 \le i < 302293`$ and $`0 \le j < 120`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Simulated annealing with targeted moves and a movable involution in the one-involution Cayley family on $`\mathbb{Z}_{302293} \rtimes_{7818} \mathbb{Z}_{120}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (235,788,540 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg13/dia8` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg13/dia8
