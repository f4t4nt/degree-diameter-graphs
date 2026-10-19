# Degree 10, Diameter 8

A graph with 4,522,590 vertices, maximum degree 10 and diameter 8; the Moore bound is 53,808,401.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{50251} \rtimes_{3} \mathbb{Z}_{90} = \langle x, y \mid x^{50251} = y^{90} = 1,\ y x y^{-1} = x^{3} \rangle`$, of order 4,522,590 (with $`50251 = 31 \cdot 1621`$), and connection set $`S = \{x^{12329} y^{24}, x^{3077} y^{35}, x^{39374} y^{34}, x^{25424} y^{39}, x^{28529} y^{26}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`S`$ contains no involution ($`G`$ has 31 involutions), so the degree is $`2 \cdot 5 = 10`$.

Vertex $`i + 50251 j`$ is $`x^i y^j`$ for $`0 \le i < 50251`$ and $`0 \le j < 90`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Even-$`M`$ group with a tiny involution class chosen by a random-set softness screen of a window of group orders, coincidence-free level pattern, simulated annealing with targeted moves on $`\mathbb{Z}_{50251} \rtimes_3 \mathbb{Z}_{90}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (22,612,950 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg10/dia8` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg10/dia8
