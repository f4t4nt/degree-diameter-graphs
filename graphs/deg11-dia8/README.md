# Degree 11, Diameter 8

A graph with 9,739,368 vertices, maximum degree 11 and diameter 8. The degree-diameter table listed 9,261,072 vertices for this cell on October 14, 2026. The Moore bound for this cell is 122,222,222.

## Construction

The Cayley graph of G = Z_135269 x|_8594 Z_72 = <x, y | x^135269 = y^72 = 1, y x y^-1 = x^8594>, of order 9,739,368 (135269 = 17 x 73 x 109), with connection set S = {x^34553 y^36, x^125066 y^35, x^6101 y^58, x^46480 y^6, x^47196 y^26, x^113742 y^19}: every g in G is joined to g s and g s^-1 for s in S.

x^34553 y^36 is one of the 1,241 involutions in G, so the degree is 1 + 2 x 5 = 11. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 135269 j is the group element x^i y^j (0 <= i < 135269, 0 <= j < 72), so vertex 0 is the identity and the 72 runs of 135,269 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 72).

## How it was found

Simulated annealing with targeted moves and a movable involution in the one-involution Cayley family on Z_135269 x|_8594 Z_72, group windows re-screened upward. The search code is in [search/cayley_metacyclic_anneal](../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (53,566,524 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg11-dia8` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg11-dia8
