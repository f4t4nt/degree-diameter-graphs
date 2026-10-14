# Degree 10, Diameter 9

A graph with 34,182,150 vertices, maximum degree 10 and diameter 9. The degree-diameter table listed 31,090,100 vertices for this cell on October 14, 2026. The Moore bound for this cell is 484,275,611.

## Construction

The Cayley graph of G = Z_227881 x|_1562 Z_150 = <x, y | x^227881 = y^150 = 1, y x y^-1 = x^1562>, of order 34,182,150 (227881 = 31 x 7351), with connection set S = {x^123970 y^23, x^17860 y^138, x^124293 y^68, x^220620 y^19, x^148855 y^3}: every g in G is joined to g s and g s^-1 for s in S.

S contains no involution (G has 31 involutions), so the degree is 2 x 5 = 10. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 227881 j is the group element x^i y^j (0 <= i < 227881, 0 <= j < 150), so vertex 0 is the identity and the 150 runs of 227,881 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 150).

## How it was found

Even-M group with a tiny involution class, simulated annealing with targeted moves from a random start on Z_227881 x|_1562 Z_150.

## Files

- No edge list is stored (170,910,750 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg10-dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg10-dia9
