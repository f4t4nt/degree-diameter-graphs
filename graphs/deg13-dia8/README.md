# Degree 13, Diameter 8

A graph with 36,275,160 vertices, maximum degree 13 and diameter 8. The degree-diameter table listed 35,140,270 vertices for this cell on October 14, 2026. The Moore bound for this cell is 508,160,186.

## Construction

The Cayley graph of G = Z_302293 x|_7818 Z_120 = <x, y | x^302293 = y^120 = 1, y x y^-1 = x^7818>, of order 36,275,160 (302293 = 41 x 73 x 101), with connection set S = {x^164125 y^60, x^129522 y^3, x^47580 y^27, x^218857 y^25, x^34159 y^56, x^366 y^18, x^71547 y^44}: every g in G is joined to g s and g s^-1 for s in S.

x^164125 y^60 is one of the 2,993 involutions in G, so the degree is 1 + 2 x 6 = 13. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 302293 j is the group element x^i y^j (0 <= i < 302293, 0 <= j < 120), so vertex 0 is the identity and the 120 runs of 302,293 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 120).

## How it was found

Simulated annealing with targeted moves and a movable involution in the one-involution Cayley family on Z_302293 x|_7818 Z_120.

## Files

- No edge list is stored (235,788,540 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg13-dia8` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg13-dia8
