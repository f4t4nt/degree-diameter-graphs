# Degree 6, Diameter 8

A graph with 80,050 vertices, maximum degree 6 and diameter 8. The Moore bound for degree 6 and diameter 8 is 585,937.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{3202} \rtimes_{1419} \mathbb{Z}_{25} = \langle x, y \mid x^{3202} = y^{25} = 1,\ y x y^{-1} = x^{1419} \rangle`$, of order 80,050 (with $`3202 = 2 \cdot 1601`$), and connection set $`S = \{x y, x y^{5}, x^{350} y^{11}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`S`$ contains no involution ($`G`$ has one involution), so the degree is $`2 \cdot 3 = 6`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 3202 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 3202`$, $`0 \le j < 25`$), so vertex 0 is the identity and the 25 runs of 3,202 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 25`$ (see [the numbering figure](../../cayley-numbering.svg)).

## How it was found

Fully symmetry-reduced complete sweep of the 3-pair Cayley graphs of the one-involution group $`\mathbb{Z}_{3202} \rtimes_{1419} \mathbb{Z}_{25}`$ (the unique diameter-8 orbit), the group chosen by the rule of few involutions and large complement order $`M`$. The search code is in [search/cayley_exhaustive](../../../search/cayley_exhaustive).

## Files

- `deg6-dia8-n80050.edges.txt`: the 240,150 edges, one `u v` per line in the numbering above, grouped by generator under `#` headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter; the edge list is also checked to be exactly the construction in `construction.json`.

    python3 tools/ddgraph.py check graphs/deg6/dia8
