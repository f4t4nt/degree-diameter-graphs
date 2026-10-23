# Degree 9, Diameter 7

A graph with 309,510 vertices, maximum degree 9 and diameter 7; the Moore bound is 2,696,338.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{3439} \rtimes_{366} \mathbb{Z}_{90} = \langle x, y \mid x^{3439} = y^{90} = 1,\ y x y^{-1} = x^{366} \rangle`$, of order 309,510 (with $`3439 = 19 \cdot 181`$), and connection set $`S = \{x^{1140} y^{45}, x^{490} y^{67}, x^{1579} y^{24}, x^{2818} y^{5}, x^{539} y^{32}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`x^{1140} y^{45}`$ is one of the 181 involutions of $`G`$, so the degree is $`1 + 2 \cdot 4 = 9`$.

Vertex $`i + 3439 j`$ is $`x^i y^j`$ for $`0 \le i < 3439`$ and $`0 \le j < 90`$, as in [the numbering figure](../../cayley-numbering.svg).

## Structure

The automorphism group is $`G`$ acting by left multiplication: the graph is a graphical regular representation of $`G`$.

## How it was found

Group chosen by a random-set softness screen of a window of group orders, simulated annealing with a fixed level pattern, targeted moves and a movable involution in the one-involution odd-degree Cayley family on $`\mathbb{Z}_{3439} \rtimes_{366} \mathbb{Z}_{90}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- `deg9-dia7-n309510.edges.txt.gz`: the 1,392,795 edges as `u v` lines, under a `#` heading per generator.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The edge list is checked against `construction.json`. The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg9/dia7
