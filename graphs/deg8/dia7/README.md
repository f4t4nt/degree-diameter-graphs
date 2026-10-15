# Degree 8, Diameter 7

A graph with 138,050 vertices, maximum degree 8 and diameter 7; the Moore bound is 1,098,057.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{2761} \rtimes_{63} \mathbb{Z}_{50} = \langle x, y \mid x^{2761} = y^{50} = 1,\ y x y^{-1} = x^{63} \rangle`$, of order 138,050 (with $`2761 = 11 \cdot 251`$), and connection set $`S = \{x^{290} y^{14}, x^{971} y^{43}, x^{457} y^{6}, x^{1728} y^{16}\}`$: each $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for each $`s \in S`$.

$`S`$ contains no involution ($`G`$ has 11 involutions), so the degree is $`2 \cdot 4 = 8`$.

Vertex $`i + 2761 j`$ is $`x^i y^j`$ for $`0 \le i < 2761`$ and $`0 \le j < 50`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Even-$`M`$ metacyclic group with a tiny involution class chosen by a random-set screen, simulated annealing followed by exact single-slot completion; $`\mathrm{Cay}(\mathbb{Z}_{2761} \rtimes_{63} \mathbb{Z}_{50})`$ with four inverse pairs. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- `deg8-dia7-n138050.edges.txt.gz`: the 552,200 edges as `u v` lines, under a `#` heading per generator.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The edge list is checked against `construction.json`. The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg8/dia7
