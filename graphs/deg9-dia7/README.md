# Degree 9, Diameter 7

A graph with 302,082 vertices, maximum degree 9 and diameter 7. The degree-diameter table listed 300,518 vertices for this cell on October 14, 2026. The Moore bound for this cell is 2,696,338.

## Construction

The Cayley graph of G = Z_4577 x|_174 Z_66 = <x, y | x^4577 = y^66 = 1, y x y^-1 = x^174>, of order 302,082 (4577 = 23 x 199), with connection set S = {y^33, x^804 y, x^3376 y^6, x^2519 y^51, x^2123 y^23}: every g in G is joined to g s and g s^-1 for s in S.

y^33 is one of the 199 involutions in G, so the degree is 1 + 2 x 4 = 9. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 4577 j is the group element x^i y^j (0 <= i < 4577, 0 <= j < 66), so vertex 0 is the identity and the 66 runs of 4,577 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 66).

## How it was found

Simulated annealing with targeted moves in the one-involution odd-degree Cayley family on Z_4577 x|_174 Z_66.

## Files

- `deg9-dia7-n302082.edges.txt.gz`: the 1,359,369 edges (gzipped), one "u v" per line in the numbering above, grouped by generator under "#" headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The edge list is checked to be exactly the construction in `construction.json`. The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg9-dia7
