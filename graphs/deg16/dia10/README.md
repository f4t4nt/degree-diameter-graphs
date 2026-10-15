# Degree 16, Diameter 10

A graph with 29,899,503,497 vertices, maximum degree 16 and diameter 10; the Moore bound is 659,029,017,857.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{24170981} \rtimes_{16573517} \mathbb{Z}_{1237} = \langle x, y \mid x^{24170981} = y^{1237} = 1,\ y x y^{-1} = x^{16573517} \rangle`$, of order 29,899,503,497 (with 24170981 prime), and connection set $`S = \{x^{10866025} y^{532}, x^{21248118} y^{539}, x^{13248079} y^{918}, x^{21841737} y^{1181}, x^{1620224} y^{1059}, x^{12379694} y^{1095}, x^{17981217} y^{470}, x^{3158481} y^{1021}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`G`$ has odd order, so it has no involutions, and the degree is $`2 \cdot 8 = 16`$.

Vertex $`i + 24170981 j`$ is $`x^i y^j`$ for $`0 \le i < 24170981`$ and $`0 \le j < 1237`$, as in [the numbering figure](../../cayley-numbering.svg).

## How it was found

Odd-order group with prime $`N`$ and $`M`$ and $`M \approx 0.40\, n^{1/3}`$, two-plane bitmap simulated annealing with targeted moves on $`\mathbb{Z}_{24170981} \rtimes_{16573517} \mathbb{Z}_{1237}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (239,196,027,976 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg16/dia10` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search, which `check` skips unless asked:

    python3 tools/ddgraph.py check graphs/deg16/dia10 --max-vertices 29899503497
