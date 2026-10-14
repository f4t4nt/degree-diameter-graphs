# Degree 12, Diameter 7

A graph with 2,114,927 vertices, maximum degree 12 and diameter 7. The degree-diameter table listed 2,031,500 vertices for this cell on October 14, 2026. The Moore bound for this cell is 23,384,605.

## Construction

The Cayley graph of G = Z_19403 x|_6308 Z_109 = <x, y | x^19403 = y^109 = 1, y x y^-1 = x^6308>, of order 2,114,927 (19403 is prime), with connection set S = {x^5982 y^44, x^1396 y^55, x^16241 y^66, x^10965 y^40, x^12766 y^46, x^5379 y^103}: every g in G is joined to g s and g s^-1 for s in S.

G has odd order, so it has no involutions, and the degree is 2 x 6 = 12. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 19403 j is the group element x^i y^j (0 <= i < 19403, 0 <= j < 109), so vertex 0 is the identity and the 109 runs of 19,403 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 109).

## How it was found

Exact single-slot completion sweep of a frozen defect-2 annealing state of 6-pair Cayley graphs on the odd-order group Z_19403 x|_6308 Z_109.

## Files

- No edge list is stored (12,689,562 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg12-dia7` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg12-dia7
