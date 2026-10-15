# Degree 13, Diameter 7

A graph with 3,691,038 vertices, maximum degree 13 and diameter 7. The degree-diameter table listed 3,529,240 vertices for this cell on October 14, 2026. The Moore bound for this cell is 42,346,682.

## Construction

The Cayley graph of G = Z_47321 x|_3981 Z_78 = <x, y | x^47321 = y^78 = 1, y x y^-1 = x^3981>, of order 3,691,038 (47321 = 79 x 599), with connection set S = {x^46136 y^39, x^28946 y^31, x^23110 y^28, x^30322 y^32, x^27810 y^11, x^25477 y^37, x^46989 y^12}: every g in G is joined to g s and g s^-1 for s in S.

x^46136 y^39 is one of the 599 involutions in G, so the degree is 1 + 2 x 6 = 13. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 47321 j is the group element x^i y^j (0 <= i < 47321, 0 <= j < 78), so vertex 0 is the identity and the 78 runs of 47,321 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 78).

## How it was found

Simulated annealing with targeted moves and a movable involution in the one-involution Cayley family on Z_47321 x|_3981 Z_78, the group chosen by the lower-quartile softness screen. The search code is in [search/cayley_metacyclic_anneal](../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (23,991,747 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg13-dia7` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg13-dia7
