# Degree 10, Diameter 8

A graph with 4,507,500 vertices, maximum degree 10 and diameter 8. The degree-diameter table listed 4,387,768 vertices for this cell on October 14, 2026. The Moore bound for this cell is 53,808,401.

## Construction

The Cayley graph of G = Z_45075 x|_2 Z_100 = <x, y | x^45075 = y^100 = 1, y x y^-1 = x^2>, of order 4,507,500 (45075 = 3 x 5^2 x 601), with connection set S = {x^35808 y^14, x^21603 y, x^8790 y^5, x^43661 y^78, x^13901 y^60}: every g in G is joined to g s and g s^-1 for s in S.

S contains no involution (G has 25 involutions), so the degree is 2 x 5 = 10. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 45075 j is the group element x^i y^j (0 <= i < 45075, 0 <= j < 100), so vertex 0 is the identity and the 100 runs of 45,075 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 100).

## How it was found

Even-M group with a tiny involution class, coincidence-free level pattern, simulated annealing with targeted moves on Z_45075 x|_2 Z_100. The search code is in [search/cayley_metacyclic_anneal](../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (22,537,500 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg10-dia8` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg10-dia8
