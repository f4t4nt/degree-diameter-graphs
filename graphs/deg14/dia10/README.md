# Degree 14, Diameter 10

A graph with 8,052,187,457 vertices, maximum degree 14 and diameter 10. The Moore bound for degree 14 and diameter 10 is 160,834,907,157.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{10837399} \rtimes_{3457} \mathbb{Z}_{743} = \langle x, y \mid x^{10837399} = y^{743} = 1,\ y x y^{-1} = x^{3457} \rangle`$, of order 8,052,187,457 (with 10837399 prime), and connection set $`S = \{x^{7548685} y^{659}, x^{870716} y^{689}, x^{375693} y^{280}, x^{7060489} y^{88}, x^{4756017} y^{500}, x^{6242423} y^{171}, x^{6945981} y^{558}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`G`$ has odd order, so it has no involutions, and the degree is $`2 \cdot 7 = 14`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 10837399 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 10837399`$, $`0 \le j < 743`$), so vertex 0 is the identity and the 743 runs of 10,837,399 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 743`$.

## How it was found

Odd-order group with $`M \approx 0.37\, n^{1/3}`$, bitmap simulated annealing with targeted moves on $`\mathbb{Z}_{10837399} \rtimes_{3457} \mathbb{Z}_{743}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (56,365,312,199 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg14/dia10` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search. It needs a few bits of memory per vertex, so `check` skips it unless asked:

    python3 tools/ddgraph.py check graphs/deg14/dia10 --max-vertices 8052187457
