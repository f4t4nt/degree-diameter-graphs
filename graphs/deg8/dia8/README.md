# Degree 8, Diameter 8

A graph with 782,838 vertices, maximum degree 8 and diameter 8; the Moore bound is 7,686,401.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{14497} \rtimes_{157} \mathbb{Z}_{54} = \langle x, y \mid x^{14497} = y^{54} = 1,\ y x y^{-1} = x^{157} \rangle`$, of order 782,838 (with $`14497 = 7 \cdot 19 \cdot 109`$), and connection set $`S = \{x^{12598} y, x^{12925} y^{11}, x^{9961} y^{50}, x^{8964} y^{24}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`S`$ contains no involution ($`G`$ has 7 involutions), so the degree is $`2 \cdot 4 = 8`$.

Vertex $`i + 14497 j`$ is $`x^i y^j`$ for $`0 \le i < 14497`$ and $`0 \le j < 54`$, as in [the numbering figure](../../cayley-numbering.svg).

## Structure

The automorphism group is $`G`$ acting by left multiplication: the graph is a graphical regular representation of $`G`$.

## How it was found

Even-$`M`$ group with a tiny involution class chosen by a random-set softness screen of a window of group orders, coincidence-free level pattern, simulated annealing with targeted moves on $`\mathbb{Z}_{14497} \rtimes_{157} \mathbb{Z}_{54}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- `deg8-dia8-n782838.edges.txt.gz`: the 3,131,352 edges as `u v` lines, under a `#` heading per generator.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The edge list is checked against `construction.json`. The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg8/dia8
