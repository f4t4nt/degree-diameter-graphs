# Degree 11, Diameter 7

A graph with 1,194,278 vertices, maximum degree 11 and diameter 7. The degree-diameter table listed 1,153,512 vertices for this cell on October 14, 2026. The Moore bound for this cell is 12,222,222.

## Construction

The Cayley graph of G = Z_20591 x|_2507 Z_58 = <x, y | x^20591 = y^58 = 1, y x y^-1 = x^2507>, of order 1,194,278 (20591 = 59 x 349), with connection set S = {x^19883 y^29, x^11052 y^10, x^2968 y^3, x^1044 y^2, x^3615 y^43, x^13568 y^37}: every g in G is joined to g s and g s^-1 for s in S.

x^19883 y^29 is one of the 349 involutions in G, so the degree is 1 + 2 x 5 = 11. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 20591 j is the group element x^i y^j (0 <= i < 20591, 0 <= j < 58), so vertex 0 is the identity and the 58 runs of 20,591 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 58).

## How it was found

Simulated annealing with level patterns, targeted moves and a movable involution in the one-involution Cayley family on Z_20591 x|_2507 Z_58.

## Files

- No edge list is stored (6,568,529 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg11-dia7` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg11-dia7
