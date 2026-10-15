# Degree 16, Diameter 8

A graph with 185,394,793 vertices, maximum degree 16 and diameter 8. The Moore bound for degree 16 and diameter 8 is 2,929,017,857.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{769273} \rtimes_{381264} \mathbb{Z}_{241} = \langle x, y \mid x^{769273} = y^{241} = 1,\ y x y^{-1} = x^{381264} \rangle`$, of order 185,394,793 (with 769273 prime), and connection set $`S = \{x^{703316} y^{10}, x^{718789} y^{97}, x^{15241} y^{58}, x^{669144} y^{17}, x^{442319} y^{21}, x^{102905} y^{45}, x^{735454} y^{9}, x^{553642} y^{74}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`G`$ has odd order, so it has no involutions, and the degree is $`2 \cdot 8 = 16`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 769273 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 769273`$, $`0 \le j < 241`$), so vertex 0 is the identity and the 241 runs of 769,273 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 241`$.

## How it was found

Odd-order group with prime complement M, bitmap simulated annealing with targeted moves on Z_769273 x|_381264 Z_241. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (1,483,158,344 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg16/dia8` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search. It needs a few bits of memory per vertex, so `check` skips it unless asked:

    python3 tools/ddgraph.py check graphs/deg16/dia8 --max-vertices 185394793
