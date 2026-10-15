# Degree 11, Diameter 10

A graph with 697,280,200 vertices, maximum degree 11 and diameter 10; the Moore bound is 12,222,222,222.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{3486401} \rtimes_{1834706} \mathbb{Z}_{200} = \langle x, y \mid x^{3486401} = y^{200} = 1,\ y x y^{-1} = x^{1834706} \rangle`$, of order 697,280,200 (with $`3486401 = 601 \cdot 5801`$), and connection set $`S = \{x^{397261} y^{100}, x^{2446299} y^{10}, x^{494688} y^{2}, x^{1969958} y^{19}, x^{3279112} y^{176}, x^{2535892} y^{47}\}`$: each $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for each $`s \in S`$.

$`x^{397261} y^{100}`$ is one of the 5,801 involutions of $`G`$, so the degree is $`1 + 2 \cdot 5 = 11`$.

Vertex $`i + 3486401 j`$ is $`x^i y^j`$ for $`0 \le i < 3486401`$ and $`0 \le j < 200`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Two-plane bitmap simulated annealing with targeted moves and a movable involution in the one-involution Cayley family on $`\mathbb{Z}_{3486401} \rtimes_{1834706} \mathbb{Z}_{200}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (3,835,041,100 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg11/dia10` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search, which `check` skips unless asked:

    python3 tools/ddgraph.py check graphs/deg11/dia10 --max-vertices 697280200
