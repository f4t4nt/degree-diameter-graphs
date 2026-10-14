# Degree 9, Diameter 10

A graph with 89,543,818 vertices, maximum degree 9 and diameter 10. The degree-diameter table listed 80,560,386 vertices for this cell on October 14, 2026. The Moore bound for this cell is 1,380,525,202.

## Construction

The Cayley graph of G = Z_491999 x|_477175 Z_182 = <x, y | x^491999 = y^182 = 1, y x y^-1 = x^477175>, of order 89,543,818 (491999 = 53 x 9283), with connection set S = {x^339836 y^91, x^406367 y, x^297390 y^5, x^419295 y^8, x^113563 y^118}: every g in G is joined to g s and g s^-1 for s in S.

x^339836 y^91 is one of the 9,283 involutions in G, so the degree is 1 + 2 x 4 = 9. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 491999 j is the group element x^i y^j (0 <= i < 491999, 0 <= j < 182), so vertex 0 is the identity and the 182 runs of 491,999 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 182).

## How it was found

Two-plane bitmap simulated annealing with targeted moves and a movable involution in the one-involution Cayley family on Z_491999 x|_477175 Z_182.

## Files

- No edge list is stored (402,947,181 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg9-dia10` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg9-dia10
