# Degree 13, Diameter 8

A graph with 36,275,160 vertices, maximum degree 13 and diameter 8. The Moore bound for degree 13 and diameter 8 is 508,160,186.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{302293} \rtimes_{7818} \mathbb{Z}_{120} = \langle x, y \mid x^{302293} = y^{120} = 1,\ y x y^{-1} = x^{7818} \rangle`$, of order 36,275,160 (with $`302293 = 41 \cdot 73 \cdot 101`$), and connection set $`S = \{x^{164125} y^{60}, x^{129522} y^{3}, x^{47580} y^{27}, x^{218857} y^{25}, x^{34159} y^{56}, x^{366} y^{18}, x^{71547} y^{44}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`x^{164125} y^{60}`$ is one of the 2,993 involutions of $`G`$, so the degree is $`1 + 2 \cdot 6 = 13`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 302293 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 302293`$, $`0 \le j < 120`$), so vertex 0 is the identity and the 120 runs of 302,293 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 120`$.

## How it was found

Simulated annealing with targeted moves and a movable involution in the one-involution Cayley family on Z_302293 x|_7818 Z_120. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (235,788,540 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg13/dia8` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg13/dia8
