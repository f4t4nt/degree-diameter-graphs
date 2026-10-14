# Degree 6, Diameter 8

A graph with 80,050 vertices, maximum degree 6 and diameter 8. The degree-diameter table listed 78,825 vertices for this cell on October 14, 2026. The Moore bound for this cell is 585,937.

## Construction

The Cayley graph of G = Z_3202 x|_1419 Z_25 = <x, y | x^3202 = y^25 = 1, y x y^-1 = x^1419>, of order 80,050 (3202 = 2 x 1601), with connection set S = {x y, x y^5, x^350 y^11}: every g in G is joined to g s and g s^-1 for s in S.

S contains no involution (G has one involution), so the degree is 2 x 3 = 6. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 3202 j is the group element x^i y^j (0 <= i < 3202, 0 <= j < 25), so vertex 0 is the identity and the 25 runs of 3,202 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 25).

## How it was found

Fully symmetry-reduced complete sweep of the 3-pair Cayley graphs of the one-involution group Z_3202 x|_1419 Z_25 (the unique diameter-8 orbit), the group chosen by the rule of few involutions and large complement order M.

## Files

- `deg6-dia8-n80050.edges.txt`: the 240,150 edges, one "u v" per line in the numbering above, grouped by generator under "#" headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter; the edge list is also checked to be exactly the construction in `construction.json`.

    python3 tools/ddgraph.py check graphs/deg6-dia8
