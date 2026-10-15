# Degree 7, Diameter 9

A graph with 1,332,774 vertices, maximum degree 7 and diameter 9. The Moore bound for degree 7 and diameter 9 is 14,108,774.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{24681} \rtimes_{788} \mathbb{Z}_{54} = \langle x, y \mid x^{24681} = y^{54} = 1,\ y x y^{-1} = x^{788} \rangle`$, of order 1,332,774 (with $`24681 = 3 \cdot 19 \cdot 433`$), and connection set $`S = \{x^{16435} y^{27}, x^{7192} y^{33}, x^{857} y^{2}, x^{4448} y^{16}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`x^{16435} y^{27}`$ is one of the 1,299 involutions of $`G`$, so the degree is $`1 + 2 \cdot 3 = 7`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 24681 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 24681`$, $`0 \le j < 54`$), so vertex 0 is the identity and the 54 runs of 24,681 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 54`$.

## How it was found

Simulated annealing with a fixed level pattern, targeted moves and a movable involution in the one-involution Cayley family on Z_24681 x|_788 Z_54. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (4,664,709 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg7/dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg7/dia9
