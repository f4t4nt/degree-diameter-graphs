# Degree 14, Diameter 9

A graph with 709,125,173 vertices, maximum degree 14 and diameter 9. The Moore bound for degree 14 and diameter 9 is 12,371,915,935.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{2104229} \rtimes_{864} \mathbb{Z}_{337} = \langle x, y \mid x^{2104229} = y^{337} = 1,\ y x y^{-1} = x^{864} \rangle`$, of order 709,125,173 (with 2104229 prime), and connection set $`S = \{x^{1325989} y^{152}, x^{2094980} y^{192}, x^{914967} y^{111}, x^{719227} y^{29}, x^{857786} y^{60}, x^{1407259} y^{16}, x^{1571238} y^{204}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`G`$ has odd order, so it has no involutions, and the degree is $`2 \cdot 7 = 14`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 2104229 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 2104229`$, $`0 \le j < 337`$), so vertex 0 is the identity and the 337 runs of 2,104,229 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 337`$.

## How it was found

Odd-order group with prime complement $`M`$ and $`M \approx 0.38\, n^{1/3}`$, bitmap simulated annealing with targeted moves on $`\mathbb{Z}_{2104229} \rtimes_{864} \mathbb{Z}_{337}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (4,963,876,211 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg14/dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search. It needs a few bits of memory per vertex, so `check` skips it unless asked:

    python3 tools/ddgraph.py check graphs/deg14/dia9 --max-vertices 709125173
