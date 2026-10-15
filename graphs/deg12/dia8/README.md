# Degree 12, Diameter 8

A graph with 19,144,241 vertices, maximum degree 12 and diameter 8. The Moore bound for degree 12 and diameter 8 is 257,230,657.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{90731} \rtimes_{61666} \mathbb{Z}_{211} = \langle x, y \mid x^{90731} = y^{211} = 1,\ y x y^{-1} = x^{61666} \rangle`$, of order 19,144,241 (with 90731 prime), and connection set $`S = \{x^{40068} y^{179}, x^{57435} y^{91}, x^{27235} y^{15}, x^{68856} y^{174}, x^{54460} y^{110}, x^{23436} y^{99}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`G`$ has odd order, so it has no involutions, and the degree is $`2 \cdot 6 = 12`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 90731 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 90731`$, $`0 \le j < 211`$), so vertex 0 is the identity and the 211 runs of 90,731 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 211`$.

## How it was found

Exact x-only completion sweep of a frozen defect-2 annealing state of 6-pair Cayley graphs on the odd-order group Z_90731 x|_61666 Z_211. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (114,865,446 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg12/dia8` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg12/dia8
