# Degree 11, Diameter 9

A graph with 81,961,120 vertices, maximum degree 11 and diameter 9. The Moore bound for degree 11 and diameter 9 is 1,222,222,222.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{512257} \rtimes_{706} \mathbb{Z}_{160} = \langle x, y \mid x^{512257} = y^{160} = 1,\ y x y^{-1} = x^{706} \rangle`$, of order 81,961,120 (with $`512257 = 97 \cdot 5281`$), and connection set $`S = \{x^{64408} y^{80}, x^{228344} y^{28}, x^{53239} y^{72}, x^{200192} y^{117}, x^{154943} y^{62}, x^{431197} y^{157}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`x^{64408} y^{80}`$ is one of the 5,281 involutions of $`G`$, so the degree is $`1 + 2 \cdot 5 = 11`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 512257 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 512257`$, $`0 \le j < 160`$), so vertex 0 is the identity and the 160 runs of 512,257 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 160`$ (see [the numbering figure](../../cayley-numbering.svg)).

## How it was found

Simulated annealing with targeted moves and a movable involution in the one-involution Cayley family on $`\mathbb{Z}_{512257} \rtimes_{706} \mathbb{Z}_{160}`$, the group chosen by the lower-quartile softness screen. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (450,786,160 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg11/dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg11/dia9
