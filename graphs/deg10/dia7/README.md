# Degree 10, Diameter 7

A graph with 624,690 vertices, maximum degree 10 and diameter 7; the Moore bound is 5,978,711.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{6941} \rtimes_{316} \mathbb{Z}_{90} = \langle x, y \mid x^{6941} = y^{90} = 1,\ y x y^{-1} = x^{316} \rangle`$, of order 624,690 (with $`6941 = 11 \cdot 631`$), and connection set $`S = \{x^{35} y^{35}, x^{509} y^{5}, x^{4285} y^{37}, x^{6911} y^{49}, x^{3639} y^{34}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`S`$ contains no involution ($`G`$ has 11 involutions), so the degree is $`2 \cdot 5 = 10`$.

Vertex $`i + 6941 j`$ is $`x^i y^j`$ for $`0 \le i < 6941`$ and $`0 \le j < 90`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Even-$`M`$ group with a tiny involution class chosen by a random-set softness screen of a window of group orders, coincidence-free level pattern, simulated annealing with targeted moves on $`\mathbb{Z}_{6941} \rtimes_{316} \mathbb{Z}_{90}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- `deg10-dia7-n624690.edges.txt.gz`: the 3,123,450 edges as `u v` lines, under a `#` heading per generator.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The edge list is checked against `construction.json`. The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg10/dia7
