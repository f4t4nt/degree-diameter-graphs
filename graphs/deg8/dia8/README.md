# Degree 8, Diameter 8

A graph with 802,890 vertices, maximum degree 8 and diameter 8; the Moore bound is 7,686,401.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{8921} \rtimes_{237} \mathbb{Z}_{90} = \langle x, y \mid x^{8921} = y^{90} = 1,\ y x y^{-1} = x^{237} \rangle`$, of order 802,890 (with $`8921 = 11 \cdot 811`$), and connection set $`S = \{x^{7422} y^{31}, x^{8137} y^{44}, x^{2812} y^{10}, x^{2784} y^{19}\}`$: each $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for each $`s \in S`$.

$`S`$ contains no involution ($`G`$ has 11 involutions), so the degree is $`2 \cdot 4 = 8`$.

Vertex $`i + 8921 j`$ is $`x^i y^j`$ for $`0 \le i < 8921`$ and $`0 \le j < 90`$, as in [the numbering figure](../../cayley-numbering.svg).

## Structure

The automorphism group is $`G`$ acting by left multiplication: the graph is a graphical regular representation of $`G`$.

## How it was found

Even-$`M`$ group with a tiny involution class chosen by a random-set softness screen of a window of group orders, coincidence-free level pattern, simulated annealing with targeted moves on $`\mathbb{Z}_{8921} \rtimes_{237} \mathbb{Z}_{90}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- `deg8-dia8-n802890.edges.txt.gz`: the 3,211,560 edges as `u v` lines, under a `#` heading per generator.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The edge list is checked against `construction.json`. The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg8/dia8
