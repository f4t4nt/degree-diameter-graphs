# Degree 12, Diameter 10

A graph with 1,692,658,087 vertices, maximum degree 12 and diameter 10; the Moore bound is 31,124,909,521.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{3655849} \rtimes_{5568} \mathbb{Z}_{463} = \langle x, y \mid x^{3655849} = y^{463} = 1,\ y x y^{-1} = x^{5568} \rangle`$, of order 1,692,658,087 (with 3655849 prime), and connection set $`S = \{x^{514157} y^{21}, x^{754332} y^{50}, x^{366326} y^{165}, x^{737423} y^{84}, x^{3525169} y^{178}, x^{3579433} y^{276}\}`$: each $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for each $`s \in S`$.

$`G`$ has odd order, so it has no involutions, and the degree is $`2 \cdot 6 = 12`$.

Vertex $`i + 3655849 j`$ is $`x^i y^j`$ for $`0 \le i < 3655849`$ and $`0 \le j < 463`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Odd-order group with small complement $`M`$, bitmap simulated annealing with targeted moves on $`\mathbb{Z}_{3655849} \rtimes_{5568} \mathbb{Z}_{463}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (10,155,948,522 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg12/dia10` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search, which `check` skips unless asked:

    python3 tools/ddgraph.py check graphs/deg12/dia10 --max-vertices 1692658087
