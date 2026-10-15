# Degree 10, Diameter 8

A graph with 4,507,500 vertices, maximum degree 10 and diameter 8. The Moore bound for degree 10 and diameter 8 is 53,808,401.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{45075} \rtimes_{2} \mathbb{Z}_{100} = \langle x, y \mid x^{45075} = y^{100} = 1,\ y x y^{-1} = x^{2} \rangle`$, of order 4,507,500 (with $`45075 = 3 \cdot 5^{2} \cdot 601`$), and connection set $`S = \{x^{35808} y^{14}, x^{21603} y, x^{8790} y^{5}, x^{43661} y^{78}, x^{13901} y^{60}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`S`$ contains no involution ($`G`$ has 25 involutions), so the degree is $`2 \cdot 5 = 10`$. The graph is vertex-transitive: $`G`$ acts on it by left multiplication.

Vertex $`i + 45075 j`$ is the group element $`x^i y^j`$ ($`0 \le i < 45075`$, $`0 \le j < 100`$), so vertex 0 is the identity and the 100 runs of 45,075 consecutive numbers are the cosets of $`\langle x \rangle`$; a generator $`x^a y^b`$ joins run $`j`$ to run $`j + b \bmod 100`$ (see [the numbering figure](../../cayley-numbering.svg)).

## How it was found

Even-$`M`$ group with a tiny involution class, coincidence-free level pattern, simulated annealing with targeted moves on $`\mathbb{Z}_{45075} \rtimes_2 \mathbb{Z}_{100}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (22,537,500 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg10/dia8` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg10/dia8
