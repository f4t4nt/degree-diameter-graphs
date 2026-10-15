# Degree 10, Diameter 7

A graph with 616,915 vertices, maximum degree 10 and diameter 7. The Moore bound for degree 10 and diameter 7 is 5,978,711.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{9491} \rtimes_{86} \mathbb{Z}_{65} = \langle x, y \mid x^{9491} = y^{65} = 1,\ y x y^{-1} = x^{86} \rangle`$, of order 616,915 (with 9491 prime), and connection set $`S = \{x^{7392} y^{64}, x^{8053} y^{11}, x^{3256} y^{30}, x^{7302} y^{41}, x^{3430} y^{9}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`G`$ has odd order, so it has no involutions, and the degree is $`2 \cdot 5 = 10`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 9491 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 9491`$, $`0 \le j < 65`$), so vertex 0 is the identity and the 65 runs of 9,491 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 65`$ (see [the numbering figure](../../cayley-numbering.svg)).

## How it was found

Odd-order group picked by a random-set screen, simulated annealing with exact single-slot completion; $`\mathrm{Cay}(\mathbb{Z}_{9491} \rtimes_{86} \mathbb{Z}_{65})`$ with five inverse pairs. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- `deg10-dia7-n616915.edges.txt.gz`: the 3,084,575 edges (gzipped), one `u v` per line in the numbering above, grouped by generator under `#` headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The edge list is checked to be exactly the construction in `construction.json`. The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg10/dia7
