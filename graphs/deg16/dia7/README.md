# Degree 16, Diameter 7

A graph with 14,995,367 vertices, maximum degree 16 and diameter 7; the Moore bound is 195,267,857.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{83773} \rtimes_{72211} \mathbb{Z}_{179} = \langle x, y \mid x^{83773} = y^{179} = 1,\ y x y^{-1} = x^{72211} \rangle`$, of order 14,995,367 (with 83773 prime), and connection set $`S = \{x^{75326} y^{15}, x^{13594} y^{16}, x^{26869} y^{58}, x^{58304} y^{33}, x^{79357} y^{12}, x^{79527} y^{69}, x^{57754} y^{71}, x^{29027} y^{7}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`G`$ has odd order, so it has no involutions, and the degree is $`2 \cdot 8 = 16`$.

Vertex $`i + 83773 j`$ is $`x^i y^j`$ for $`0 \le i < 83773`$ and $`0 \le j < 179`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Coincidence-free level pattern from the pattern optimiser, then simulated annealing with targeted moves on the odd-order group $`\mathbb{Z}_{83773} \rtimes_{72211} \mathbb{Z}_{179}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (119,962,936 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg16/dia7` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg16/dia7
