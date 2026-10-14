# Degree 15, Diameter 7

A graph with 9,758,988 vertices, maximum degree 15 and diameter 7. The degree-diameter table listed 9,372,780 vertices for this cell on October 14, 2026. The Moore bound for this cell is 121,630,966.

## Construction

The Cayley graph of G = Z_90361 x|_1443 Z_108 = <x, y | x^90361 = y^108 = 1, y x y^-1 = x^1443>, of order 9,758,988 (90361 = 109 x 829), with connection set S = {x^73139 y^54, x^62125 y^29, x^25137 y^47, x^26271 y^26, x^24758 y^37, x^18988 y^31, x^31 y^11, x^14702 y^70}: every g in G is joined to g s and g s^-1 for s in S.

x^73139 y^54 is one of the 829 involutions in G, so the degree is 1 + 2 x 7 = 15. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 90361 j is the group element x^i y^j (0 <= i < 90361, 0 <= j < 108), so vertex 0 is the identity and the 108 runs of 90,361 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 108).

## How it was found

Simulated annealing with targeted moves and a movable involution in the one-involution Cayley family on Z_90361 x|_1443 Z_108, the group chosen by the lower-quartile softness screen.

## Files

- No edge list is stored (73,192,410 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg15-dia7` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg15-dia7
