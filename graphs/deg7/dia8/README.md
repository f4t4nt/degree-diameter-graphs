# Degree 7, Diameter 8

A graph with 264,440 vertices, maximum degree 7 and diameter 8. The degree-diameter table listed 255,938 vertices for this cell on October 14, 2026. The Moore bound for this cell is 2,351,462.

## Construction

The Cayley graph of G = Z_6611 x|_193 Z_40 = <x, y | x^6611 = y^40 = 1, y x y^-1 = x^193>, of order 264,440 (6611 = 11 x 601), with connection set S = {y^20, x y^4, x^4755 y^9, x^2638 y^6}: every g in G is joined to g s and g s^-1 for s in S.

y^20 is one of the 601 involutions in G, so the degree is 1 + 2 x 3 = 7. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 6611 j is the group element x^i y^j (0 <= i < 6611, 0 <= j < 40), so vertex 0 is the identity and the 40 runs of 6,611 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 40).

## How it was found

Complete mode-F sweeps of the unit-x one-involution Cayley family per (group, level pattern) with the automorphism quotient; Cay(Z_6611 x|_193 Z_40) with the involution y^20 and three inverse pairs. The search code is in [search/cayley_exhaustive](../../search/cayley_exhaustive).

## Files

- `deg7-dia8-n264440.edges.txt.gz`: the 925,540 edges (gzipped), one "u v" per line in the numbering above, grouped by generator under "#" headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The edge list is checked to be exactly the construction in `construction.json`. The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg7-dia8
