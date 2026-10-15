# Degree 14, Diameter 8

A graph with 64,298,753 vertices, maximum degree 14 and diameter 8. The degree-diameter table listed 58,673,788 vertices for this cell on October 14, 2026. The Moore bound for this cell is 951,685,841.

## Construction

The Cayley graph of G = Z_185299 x|_47110 Z_347 = <x, y | x^185299 = y^347 = 1, y x y^-1 = x^47110>, of order 64,298,753 (185299 is prime), with connection set S = {x^101500 y^2, x^120355 y^7, x^99761 y^10, x^113037 y^111, x^146748 y^25, x^71617 y^144, x^184902 y^291}: every g in G is joined to g s and g s^-1 for s in S.

G has odd order, so it has no involutions, and the degree is 2 x 7 = 14. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 185299 j is the group element x^i y^j (0 <= i < 185299, 0 <= j < 347), so vertex 0 is the identity and the 347 runs of 185,299 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 347).

## How it was found

Odd-order group, coincidence-free level pattern with the lowest yardstick among the pattern-optimiser restarts, simulated annealing with targeted moves on Z_185299 x|_47110 Z_347. The search code is in [search/cayley_metacyclic_anneal](../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (450,091,271 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg14-dia8` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg14-dia8
