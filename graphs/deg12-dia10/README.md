# Degree 12, Diameter 10

A graph with 1,692,658,087 vertices, maximum degree 12 and diameter 10. The degree-diameter table listed 1,545,004,522 vertices for this cell on October 14, 2026. The Moore bound for this cell is 31,124,909,521.

## Construction

The Cayley graph of G = Z_3655849 x|_5568 Z_463 = <x, y | x^3655849 = y^463 = 1, y x y^-1 = x^5568>, of order 1,692,658,087 (3655849 is prime), with connection set S = {x^514157 y^21, x^754332 y^50, x^366326 y^165, x^737423 y^84, x^3525169 y^178, x^3579433 y^276}: every g in G is joined to g s and g s^-1 for s in S.

G has odd order, so it has no involutions, and the degree is 2 x 6 = 12. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 3655849 j is the group element x^i y^j (0 <= i < 3655849, 0 <= j < 463), so vertex 0 is the identity and the 463 runs of 3,655,849 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 463).

## How it was found

Odd-order group with small complement M, bitmap simulated annealing with targeted moves on Z_3655849 x|_5568 Z_463.

## Files

- No edge list is stored (10,155,948,522 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg12-dia10` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search. It needs a few bits of memory per vertex, so `check` skips it unless asked:

    python3 tools/ddgraph.py check graphs/deg12-dia10 --max-vertices 1692658087
