# Degree 8, Diameter 9

A graph with 4,601,361 vertices, maximum degree 8 and diameter 9; the Moore bound is 53,804,809.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{49477} \rtimes_{19070} \mathbb{Z}_{93} = \langle x, y \mid x^{49477} = y^{93} = 1,\ y x y^{-1} = x^{19070} \rangle`$, of order 4,601,361 (with 49477 prime), and connection set $`S = \{x^{24458} y^{25}, x^{25770} y^{85}, x^{7723} y^{21}, x^{40638} y^{79}\}`$: each $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for each $`s \in S`$.

$`G`$ has odd order, so it has no involutions, and the degree is $`2 \cdot 4 = 8`$.

Vertex $`i + 49477 j`$ is $`x^i y^j`$ for $`0 \le i < 49477`$ and $`0 \le j < 93`$, as in [the numbering figure](../../cayley-numbering.svg).

## Structure

The automorphism group is $`G`$ acting by left multiplication: the graph is a graphical regular representation of $`G`$.

## How it was found

Odd-order group chosen by a random-set softness screen of a window of group orders, coincidence-free level pattern, simulated annealing with targeted moves on $`\mathbb{Z}_{49477} \rtimes_{19070} \mathbb{Z}_{93}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (18,405,444 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg8/dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg8/dia9
