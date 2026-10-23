# Degree 8, Diameter 9

A graph with 4,516,807 vertices, maximum degree 8 and diameter 9; the Moore bound is 53,804,809.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{63617} \rtimes_{317} \mathbb{Z}_{71} = \langle x, y \mid x^{63617} = y^{71} = 1,\ y x y^{-1} = x^{317} \rangle`$, of order 4,516,807 (with 63617 prime), and connection set $`S = \{x^{56818} y^{30}, x^{57247} y^{28}, x^{27552} y^{16}, x^{8832} y^{34}\}`$: each $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for each $`s \in S`$.

$`G`$ has odd order, so it has no involutions, and the degree is $`2 \cdot 4 = 8`$.

Vertex $`i + 63617 j`$ is $`x^i y^j`$ for $`0 \le i < 63617`$ and $`0 \le j < 71`$, as in [the numbering figure](../../cayley-numbering.svg).

## Structure

The automorphism group is $`G`$ acting by left multiplication: the graph is a graphical regular representation of $`G`$.

## How it was found

Odd-order group chosen by a random-set softness screen of a window of group orders, coincidence-free level pattern, simulated annealing with targeted moves followed by exact single-slot completion on $`\mathbb{Z}_{63617} \rtimes_{317} \mathbb{Z}_{71}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (18,067,228 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg8/dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg8/dia9
