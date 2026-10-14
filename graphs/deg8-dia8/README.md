# Degree 8, Diameter 8

A graph with 779,803 vertices, maximum degree 8 and diameter 8. The degree-diameter table listed 744,942 vertices for this cell on October 14, 2026. The Moore bound for this cell is 7,686,401.

## Construction

The Cayley graph of G = Z_13217 x|_8875 Z_59 = <x, y | x^13217 = y^59 = 1, y x y^-1 = x^8875>, of order 779,803 (13217 is prime), with connection set S = {x^761 y, x^3375 y^12, x^11838 y^37, x^5435 y^28}: every g in G is joined to g s and g s^-1 for s in S.

G has odd order, so it has no involutions, and the degree is 2 x 4 = 8. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 13217 j is the group element x^i y^j (0 <= i < 13217, 0 <= j < 59), so vertex 0 is the identity and the 59 runs of 13,217 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 59).

## How it was found

Coincidence-free level pattern from the pattern optimiser, then simulated annealing with targeted moves on the odd-order group Z_13217 x|_8875 Z_59.

## Files

- `deg8-dia8-n779803.edges.txt.gz`: the 3,119,212 edges (gzipped), one "u v" per line in the numbering above, grouped by generator under "#" headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The edge list is checked to be exactly the construction in `construction.json`. The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg8-dia8
