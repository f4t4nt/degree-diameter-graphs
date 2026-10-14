# Degree 16, Diameter 7

A graph with 14,995,367 vertices, maximum degree 16 and diameter 7. The degree-diameter table listed 14,882,658 vertices for this cell on October 14, 2026. The Moore bound for this cell is 195,267,857.

## Construction

The Cayley graph of G = Z_83773 x|_72211 Z_179 = <x, y | x^83773 = y^179 = 1, y x y^-1 = x^72211>, of order 14,995,367 (83773 is prime), with connection set S = {x^75326 y^15, x^13594 y^16, x^26869 y^58, x^58304 y^33, x^79357 y^12, x^79527 y^69, x^57754 y^71, x^29027 y^7}: every g in G is joined to g s and g s^-1 for s in S.

G has odd order, so it has no involutions, and the degree is 2 x 8 = 16. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 83773 j is the group element x^i y^j (0 <= i < 83773, 0 <= j < 179), so vertex 0 is the identity and the 179 runs of 83,773 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 179).

## How it was found

Coincidence-free level pattern from the pattern optimiser, then simulated annealing with targeted moves on the odd-order group Z_83773 x|_72211 Z_179.

## Files

- No edge list is stored (119,962,936 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg16-dia7` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg16-dia7
