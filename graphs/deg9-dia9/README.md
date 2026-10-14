# Degree 9, Diameter 9

A graph with 13,219,318 vertices, maximum degree 9 and diameter 9. The degree-diameter table listed 12,441,690 vertices for this cell on October 14, 2026. The Moore bound for this cell is 172,565,650.

## Construction

The Cayley graph of G = Z_134891 x|_4746 Z_98 = <x, y | x^134891 = y^98 = 1, y x y^-1 = x^4746>, of order 13,219,318 (134891 = 43 x 3137), with connection set S = {x^48203 y^49, x^74247 y, x^51808 y^79, x^121079 y^27, x^71015 y^82}: every g in G is joined to g s and g s^-1 for s in S.

x^48203 y^49 is one of the 3,137 involutions in G, so the degree is 1 + 2 x 4 = 9. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 134891 j is the group element x^i y^j (0 <= i < 134891, 0 <= j < 98), so vertex 0 is the identity and the 98 runs of 134,891 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 98).

## How it was found

Simulated annealing with targeted moves and a movable involution in the one-involution Cayley family on Z_134891 x|_4746 Z_98, the group chosen by the lower-quartile softness screen.

## Files

- No edge list is stored (59,486,931 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg9-dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg9-dia9
