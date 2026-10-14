# Degree 9, Diameter 5

A graph with 9,020 vertices, maximum degree 9 and diameter 5. The degree-diameter table listed 8,856 vertices for this cell on October 14, 2026. The Moore bound for this cell is 42,130.

## Construction

The Cayley graph of G = Z_451 x|_46 Z_20 = <x, y | x^451 = y^20 = 1, y x y^-1 = x^46>, of order 9,020 (451 = 11 x 41), with connection set S = {y^10, x y, x^316 y^6, x^13 y^7, x^426 y^9}: every g in G is joined to g s and g s^-1 for s in S.

y^10 is one of the 41 involutions in G, so the degree is 1 + 2 x 4 = 9. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 451 j is the group element x^i y^j (0 <= i < 451, 0 <= j < 20), so vertex 0 is the identity and the 20 runs of 451 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 20).

## How it was found

Complete chunked sweep of the one-involution Cayley family {y^10, four inverse pairs} on the groups Z_451 x|_s Z_20 up to automorphisms (two inequivalent solutions on two twists s).

## Files

- `deg9-dia5-n9020.edges.txt`: the 40,590 edges, one "u v" per line in the numbering above, grouped by generator under "#" headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter; the edge list is also checked to be exactly the construction in `construction.json`.

    python3 tools/ddgraph.py check graphs/deg9-dia5
