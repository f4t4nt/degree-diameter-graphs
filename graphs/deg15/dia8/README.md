# Degree 15, Diameter 8

A graph with 110,115,980 vertices, maximum degree 15 and diameter 8. The Moore bound for degree 15 and diameter 8 is 1,702,833,526.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{423523} \rtimes_{8} \mathbb{Z}_{260} = \langle x, y \mid x^{423523} = y^{260} = 1,\ y x y^{-1} = x^{8} \rangle`$, of order 110,115,980 (with $`423523 = 53 \cdot 61 \cdot 131`$), and connection set $`S = \{x^{295274} y^{130}, x^{422256} y^{171}, x^{6132} y^{15}, x^{168289} y^{242}, x^{398805} y^{70}, x^{25253} y^{246}, x^{329528} y^{8}, x^{347373} y^{113}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`x^{295274} y^{130}`$ is one of the 3,233 involutions of $`G`$, so the degree is $`1 + 2 \cdot 7 = 15`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 423523 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 423523`$, $`0 \le j < 260`$), so vertex 0 is the identity and the 260 runs of 423,523 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 260`$ (see [the numbering figure](../../cayley-numbering.svg)).

## How it was found

Whole-window small-complement screen, simulated annealing with involution moves in the one-involution Cayley family on $`\mathbb{Z}_{423523} \rtimes_8 \mathbb{Z}_{260}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (825,869,850 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg15/dia8` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search. It needs a few bits of memory per vertex, so `check` skips it unless asked:

    python3 tools/ddgraph.py check graphs/deg15/dia8 --max-vertices 110115980
