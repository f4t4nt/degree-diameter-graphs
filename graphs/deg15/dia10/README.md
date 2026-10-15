# Degree 15, Diameter 10

A graph with 15,810,645,560 vertices, maximum degree 15 and diameter 10; the Moore bound is 333,755,371,126.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{20803481} \rtimes_{5726682} \mathbb{Z}_{760} = \langle x, y \mid x^{20803481} = y^{760} = 1,\ y x y^{-1} = x^{5726682} \rangle`$, of order 15,810,645,560 (with $`20803481 = 3041 \cdot 6841`$), and connection set $`S = \{x^{4330384} y^{380}, x^{16555534} y^{620}, x^{11877752} y^{102}, x^{20283801} y^{421}, x^{2712392} y^{488}, x^{8353947} y^{95}, x^{9367781} y^{688}, x^{5818265} y^{505}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`x^{4330384} y^{380}`$ is one of the 6,841 involutions of $`G`$, so the degree is $`1 + 2 \cdot 7 = 15`$.

Vertex $`i + 20803481 j`$ is $`x^i y^j`$ for $`0 \le i < 20803481`$ and $`0 \le j < 760`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Two-plane bitmap simulated annealing with targeted moves and a movable involution in the one-involution Cayley family on $`\mathbb{Z}_{20803481} \rtimes_{5726682} \mathbb{Z}_{760}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (118,579,841,700 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg15/dia10` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search, which `check` skips unless asked:

    python3 tools/ddgraph.py check graphs/deg15/dia10 --max-vertices 15810645560
