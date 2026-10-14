# Degree 12, Diameter 8

A graph with 19,144,241 vertices, maximum degree 12 and diameter 8. The degree-diameter table listed 18,101,302 vertices for this cell on October 14, 2026. The Moore bound for this cell is 257,230,657.

## Construction

The Cayley graph of G = Z_90731 x|_61666 Z_211 = <x, y | x^90731 = y^211 = 1, y x y^-1 = x^61666>, of order 19,144,241 (90731 is prime), with connection set S = {x^40068 y^179, x^57435 y^91, x^27235 y^15, x^68856 y^174, x^54460 y^110, x^23436 y^99}: every g in G is joined to g s and g s^-1 for s in S.

G has odd order, so it has no involutions, and the degree is 2 x 6 = 12. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 90731 j is the group element x^i y^j (0 <= i < 90731, 0 <= j < 211), so vertex 0 is the identity and the 211 runs of 90,731 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 211).

## How it was found

Exact x-only completion sweep of a frozen defect-2 annealing state of 6-pair Cayley graphs on the odd-order group Z_90731 x|_61666 Z_211.

## Files

- No edge list is stored (114,865,446 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg12-dia8` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg12-dia8
