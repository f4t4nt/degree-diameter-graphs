# Degree 10, Diameter 9

A graph with 34,182,150 vertices, maximum degree 10 and diameter 9; the Moore bound is 484,275,611.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{227881} \rtimes_{1562} \mathbb{Z}_{150} = \langle x, y \mid x^{227881} = y^{150} = 1,\ y x y^{-1} = x^{1562} \rangle`$, of order 34,182,150 (with $`227881 = 31 \cdot 7351`$), and connection set $`S = \{x^{123970} y^{23}, x^{17860} y^{138}, x^{124293} y^{68}, x^{220620} y^{19}, x^{148855} y^{3}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`S`$ contains no involution ($`G`$ has 31 involutions), so the degree is $`2 \cdot 5 = 10`$.

Vertex $`i + 227881 j`$ is $`x^i y^j`$ for $`0 \le i < 227881`$ and $`0 \le j < 150`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Even-$`M`$ group with a tiny involution class, simulated annealing with targeted moves from a random start on $`\mathbb{Z}_{227881} \rtimes_{1562} \mathbb{Z}_{150}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (170,910,750 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg10/dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg10/dia9
