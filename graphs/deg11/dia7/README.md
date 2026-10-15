# Degree 11, Diameter 7

A graph with 1,194,278 vertices, maximum degree 11 and diameter 7; the Moore bound is 12,222,222.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{20591} \rtimes_{2507} \mathbb{Z}_{58} = \langle x, y \mid x^{20591} = y^{58} = 1,\ y x y^{-1} = x^{2507} \rangle`$, of order 1,194,278 (with $`20591 = 59 \cdot 349`$), and connection set $`S = \{x^{19883} y^{29}, x^{11052} y^{10}, x^{2968} y^{3}, x^{1044} y^{2}, x^{3615} y^{43}, x^{13568} y^{37}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`x^{19883} y^{29}`$ is one of the 349 involutions of $`G`$, so the degree is $`1 + 2 \cdot 5 = 11`$.

Vertex $`i + 20591 j`$ is $`x^i y^j`$ for $`0 \le i < 20591`$ and $`0 \le j < 58`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Simulated annealing with level patterns, targeted moves and a movable involution in the one-involution Cayley family on $`\mathbb{Z}_{20591} \rtimes_{2507} \mathbb{Z}_{58}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (6,568,529 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg11/dia7` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg11/dia7
