# Degree 13, Diameter 9

A graph with 362,492,640 vertices, maximum degree 13 and diameter 9; the Moore bound is 6,097,922,234.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{1258655} \rtimes_{263} \mathbb{Z}_{288} = \langle x, y \mid x^{1258655} = y^{288} = 1,\ y x y^{-1} = x^{263} \rangle`$, of order 362,492,640 (with $`1258655 = 5 \cdot 19 \cdot 13249`$), and connection set $`S = \{x^{286140} y^{144}, x^{527543} y^{2}, x^{579807} y^{69}, x^{306882} y^{14}, x^{402861} y^{32}, x^{244866} y^{149}, x^{463227} y^{58}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`x^{286140} y^{144}`$ is one of the 13,249 involutions of $`G`$, so the degree is $`1 + 2 \cdot 6 = 13`$.

Vertex $`i + 1258655 j`$ is $`x^i y^j`$ for $`0 \le i < 1258655`$ and $`0 \le j < 288`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Small-complement group screen, bitmap simulated annealing with a movable involution in the one-involution Cayley family on $`\mathbb{Z}_{1258655} \rtimes_{263} \mathbb{Z}_{288}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (2,356,202,160 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg13/dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search, which `check` skips unless asked:

    python3 tools/ddgraph.py check graphs/deg13/dia9 --max-vertices 362492640
