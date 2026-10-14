# Degree 7, Diameter 6

A graph with 12,654 vertices, maximum degree 7 and diameter 6. The degree-diameter table listed 12,642 vertices for this cell on October 14, 2026. The Moore bound for this cell is 65,318.

## Construction

The Cayley graph of G = Z_703 x|_85 Z_18 = <x, y | x^703 = y^18 = 1, y x y^-1 = x^85>, of order 12,654 (703 = 19 x 37), with connection set S = {y^9, x y, x^24 y^5, x^289 y^8}: every g in G is joined to g s and g s^-1 for s in S.

y^9 is one of the 37 involutions in G, so the degree is 1 + 2 x 3 = 7. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 703 j is the group element x^i y^j (0 <= i < 703, 0 <= j < 18), so vertex 0 is the identity and the 18 runs of 703 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 18).

## How it was found

Complete sweep of all one-involution connection sets {y^9, three inverse pairs} on Z_703 x|_85 Z_18 up to the automorphisms fixing the involution (exactly one solution).

## Files

- `deg7-dia6-n12654.edges.txt`: the 44,289 edges, one "u v" per line in the numbering above, grouped by generator under "#" headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter; the edge list is also checked to be exactly the construction in `construction.json`.

    python3 tools/ddgraph.py check graphs/deg7-dia6
