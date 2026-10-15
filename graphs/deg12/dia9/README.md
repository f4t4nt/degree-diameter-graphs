# Degree 12, Diameter 9

A graph with 177,775,277 vertices, maximum degree 12 and diameter 9. The Moore bound for degree 12 and diameter 9 is 2,829,537,229.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{783151} \rtimes_{1808} \mathbb{Z}_{227} = \langle x, y \mid x^{783151} = y^{227} = 1,\ y x y^{-1} = x^{1808} \rangle`$, of order 177,775,277 (with 783151 prime), and connection set $`S = \{x^{262378} y^{29}, x^{52438} y^{90}, x^{51027} y^{186}, x^{163372} y^{33}, x^{6796} y^{10}, x^{575187} y^{193}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`G`$ has odd order, so it has no involutions, and the degree is $`2 \cdot 6 = 12`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 783151 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 783151`$, $`0 \le j < 227`$), so vertex 0 is the identity and the 227 runs of 783,151 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 227`$ (see [the numbering figure](../../cayley-numbering.svg)).

## How it was found

Groups chosen by the small-complement rule, bitmap simulated annealing with targeted moves on the odd-order group $`\mathbb{Z}_{783151} \rtimes_{1808} \mathbb{Z}_{227}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (1,066,651,662 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg12/dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search. It needs a few bits of memory per vertex, so `check` skips it unless asked:

    python3 tools/ddgraph.py check graphs/deg12/dia9 --max-vertices 177775277
