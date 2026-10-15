# Degree 15, Diameter 9

A graph with 1,315,107,690 vertices, maximum degree 15 and diameter 9. The Moore bound for degree 15 and diameter 9 is 23,839,669,366.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{3372071} \rtimes_{832} \mathbb{Z}_{390} = \langle x, y \mid x^{3372071} = y^{390} = 1,\ y x y^{-1} = x^{832} \rangle`$, of order 1,315,107,690 (with $`3372071 = 131 \cdot 25741`$), and connection set $`S = \{x^{959706} y^{195}, x^{810691} y, x^{2904441} y^{104}, x^{1086462} y^{18}, x^{499646} y^{4}, x^{2342996} y^{253}, x^{3240389} y^{343}, x^{399444} y^{359}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`x^{959706} y^{195}`$ is one of the 25,741 involutions of $`G`$, so the degree is $`1 + 2 \cdot 7 = 15`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 3372071 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 3372071`$, $`0 \le j < 390`$), so vertex 0 is the identity and the 390 runs of 3,372,071 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 390`$ (see [the numbering figure](../../cayley-numbering.svg)).

## How it was found

Small-complement group screen, bitmap simulated annealing with a movable involution in the one-involution Cayley family on $`\mathbb{Z}_{3372071} \rtimes_{832} \mathbb{Z}_{390}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (9,863,307,675 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg15/dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search. It needs a few bits of memory per vertex, so `check` skips it unless asked:

    python3 tools/ddgraph.py check graphs/deg15/dia9 --max-vertices 1315107690
