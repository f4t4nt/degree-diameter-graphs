# Degree 12, Diameter 9

A graph with 179,227,730 vertices, maximum degree 12 and diameter 9; the Moore bound is 2,829,537,229.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{779251} \rtimes_{1183} \mathbb{Z}_{230} = \langle x, y \mid x^{779251} = y^{230} = 1,\ y x y^{-1} = x^{1183} \rangle`$, of order 179,227,730 (with $`779251 = 11 \cdot 70841`$), and connection set $`S = \{x^{424056} y^{63}, x^{624662} y^{15}, x^{429226} y^{31}, x^{2928} y^{88}, x^{523926} y^{122}, x^{660361} y^{164}\}`$: every $`g \in G`$ is adjacent to $`g s`$ and $`g s^{-1}`$ for every $`s \in S`$.

$`S`$ contains no involution ($`G`$ has 11 involutions), so the degree is $`2 \cdot 6 = 12`$.

Vertex $`i + 779251 j`$ is $`x^i y^j`$ for $`0 \le i < 779251`$ and $`0 \le j < 230`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Even-$`M`$ group with a tiny involution class chosen by a random-set softness screen of a window of group orders, coincidence-free level pattern, simulated annealing with targeted moves on $`\mathbb{Z}_{779251} \rtimes_{1183} \mathbb{Z}_{230}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (1,075,366,380 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg12/dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search, which `check` skips unless asked:

    python3 tools/ddgraph.py check graphs/deg12/dia9 --max-vertices 179227730
