# Degree 6, Diameter 10

A graph with 1,394,189 vertices, maximum degree 6 and diameter 10. The degree-diameter table listed 1,279,794 vertices for this cell on October 14, 2026. The Moore bound for this cell is 14,648,437.

## Construction

The Cayley graph of G = Z_32423 x|_1693 Z_43 = <x, y | x^32423 = y^43 = 1, y x y^-1 = x^1693>, of order 1,394,189 (32423 is prime), with connection set S = {y^2, x y^13, x^4406 y^18}: every g in G is joined to g s and g s^-1 for s in S.

G has odd order, so it has no involutions, and the degree is 2 x 3 = 6. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 32423 j is the group element x^i y^j (0 <= i < 32423, 0 <= j < 43), so vertex 0 is the identity and the 43 runs of 32,423 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 43).

## How it was found

Complete level-pattern sweeps (every x-assignment of every level pattern) of 3-pair Cayley graphs on odd-order metacyclic groups; the graph is Cay(Z_32423 x|_1693 Z_43; y^2, x y^13, x^4406 y^18). The search code is in [search/cayley_exhaustive](../../search/cayley_exhaustive).

## Files

- No edge list is stored (4,182,567 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg6-dia10` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg6-dia10
