# Degree 9, Diameter 7

A graph with 306,852 vertices, maximum degree 9 and diameter 7; the Moore bound is 2,696,338.

## Construction

The Cayley graph of $`G = \mathbb{Z}_{3653} \rtimes_{2006} \mathbb{Z}_{84} = \langle x, y \mid x^{3653} = y^{84} = 1,\ y x y^{-1} = x^{2006} \rangle`$, of order 306,852 (with $`3653 = 13 \cdot 281`$), and connection set $`S = \{x^{2561} y^{42}, x^{957} y^{6}, x^{979} y^{41}, x^{1991} y^{62}, x^{2138} y^{17}\}`$: every $`g \in G`$ is joined to $`g s`$ and $`g s^{-1}`$ for $`s \in S`$.

$`x^{2561} y^{42}`$ is one of the 281 involutions of $`G`$, so the degree is $`1 + 2 \cdot 4 = 9`$.

Vertex $`i + 3653 j`$ is $`x^i y^j`$ for $`0 \le i < 3653`$ and $`0 \le j < 84`$, as in [the numbering figure](../../cayley-numbering.svg).

## Structure

The automorphism group is $`G`$ acting by left multiplication: the graph is a graphical regular representation of $`G`$.

## How it was found

Group chosen by a random-set softness screen of a window of group orders, simulated annealing with a fixed level pattern, targeted moves and a movable involution in the one-involution odd-degree Cayley family on $`\mathbb{Z}_{3653} \rtimes_{2006} \mathbb{Z}_{84}`$. The search code is in [search/cayley_metacyclic_anneal](../../../search/cayley_metacyclic_anneal).

## Files

- `deg9-dia7-n306852.edges.txt.gz`: the 1,380,834 edges as `u v` lines, under a `#` heading per generator.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The edge list is checked against `construction.json`. The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg9/dia7
