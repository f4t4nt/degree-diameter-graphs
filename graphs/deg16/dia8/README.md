# Degree 16, Diameter 8

A graph with 185,394,793 vertices, maximum degree 16 and diameter 8. The degree-diameter table listed 162,001,414 vertices for this cell on October 14, 2026. The Moore bound for this cell is 2,929,017,857.

## Construction

The Cayley graph of G = Z_769273 x|_381264 Z_241 = <x, y | x^769273 = y^241 = 1, y x y^-1 = x^381264>, of order 185,394,793 (769273 is prime), with connection set S = {x^703316 y^10, x^718789 y^97, x^15241 y^58, x^669144 y^17, x^442319 y^21, x^102905 y^45, x^735454 y^9, x^553642 y^74}: every g in G is joined to g s and g s^-1 for s in S.

G has odd order, so it has no involutions, and the degree is 2 x 8 = 16. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 769273 j is the group element x^i y^j (0 <= i < 769273, 0 <= j < 241), so vertex 0 is the identity and the 241 runs of 769,273 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 241).

## How it was found

Odd-order group with prime complement M, bitmap simulated annealing with targeted moves on Z_769273 x|_381264 Z_241. The search code is in [search/cayley_metacyclic_anneal](../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (1,483,158,344 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg16-dia8` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search. It needs a few bits of memory per vertex, so `check` skips it unless asked:

    python3 tools/ddgraph.py check graphs/deg16-dia8 --max-vertices 185394793
