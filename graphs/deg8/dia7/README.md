# Degree 8, Diameter 7

A graph with 139,679 vertices, maximum degree 8 and diameter 7; the Moore bound is 1,098,057.

## Construction

The Cayley graph of the group $`G = \mathbb{Z}_{6073} \rtimes_{2685} \mathbb{Z}_{23} = \langle x, y \mid x^{6073} = y^{23} = 1,\ y x y^{-1} = x^{2685} \rangle`$, of order 139,679 (with 6073 prime), and connection set $`S = \{x^{5228} y^{2}, x^{5308} y^{16}, x^{987} y^{14}, x^{5277} y^{17}\}`$: each $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for each $`s \in S`$.

$`G`$ has odd order, so it has no involutions, and the degree is $`2 \cdot 4 = 8`$.

Vertex $`i + 6073 j`$ is $`x^i y^j`$ for $`0 \le i < 6073`$ and $`0 \le j < 23`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Odd-order group chosen by a random-set softness screen of a window of group orders, the pattern optimiser's level pattern, simulated annealing with targeted moves on $`\mathbb{Z}_{6073} \rtimes_{2685} \mathbb{Z}_{23}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- `deg8-dia7-n139679.edges.txt.gz`: the 558,716 edges as `u v` lines, under a `#` heading per generator.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The edge list is checked against `construction.json`. The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg8/dia7
