# Degree 7, Diameter 10

A graph with 6,641,680 vertices, maximum degree 7 and diameter 10. The Moore bound for degree 7 and diameter 10 is 84,652,646.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{83021} \rtimes_{4400} \mathbb{Z}_{80} = \langle x, y \mid x^{83021} = y^{80} = 1,\ y x y^{-1} = x^{4400} \rangle`$, of order 6,641,680 (with $`83021 = 61 \cdot 1361`$), and connection set $`S = \{x^{74603} y^{40}, x^{4907} y^{38}, x^{21040} y^{3}, x^{63695} y^{58}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`x^{74603} y^{40}`$ is one of the 1,361 involutions of $`G`$, so the degree is $`1 + 2 \cdot 3 = 7`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 83021 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 83021`$, $`0 \le j < 80`$), so vertex 0 is the identity and the 80 runs of 83,021 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 80`$.

## How it was found

Simulated annealing in the one-involution Cayley family with involution moves that conjugate the connection set towards an unreached involution, on Z_83021 x|_4400 Z_80. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (23,245,880 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg7/dia10` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg7/dia10
