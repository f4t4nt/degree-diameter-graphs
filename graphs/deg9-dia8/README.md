# Degree 9, Diameter 8

A graph with 1,961,370 vertices, maximum degree 9 and diameter 8. The degree-diameter table listed 1,859,260 vertices for this cell on October 14, 2026. The Moore bound for this cell is 21,570,706.

## Construction

The Cayley graph of G = Z_21793 x|_34 Z_90 = <x, y | x^21793 = y^90 = 1, y x y^-1 = x^34>, of order 1,961,370 (21793 = 19 x 31 x 37), with connection set S = {x^18500 y^45, x^10329 y^8, x^6476 y^9, x^17755 y^5, x^3186 y^39}: every g in G is joined to g s and g s^-1 for s in S.

x^18500 y^45 is one of the 589 involutions in G, so the degree is 1 + 2 x 4 = 9. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 21793 j is the group element x^i y^j (0 <= i < 21793, 0 <= j < 90), so vertex 0 is the identity and the 90 runs of 21,793 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 90).

## How it was found

Simulated annealing with level patterns, targeted moves and a movable involution in the one-involution Cayley family on Z_21793 x|_34 Z_90. The search code is in [search/cayley_metacyclic_anneal](../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (8,826,165 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg9-dia8` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg9-dia8
