# Degree 13, Diameter 9

A graph with 362,492,640 vertices, maximum degree 13 and diameter 9. The degree-diameter table listed 334,485,050 vertices for this cell on October 14, 2026. The Moore bound for this cell is 6,097,922,234.

## Construction

The Cayley graph of G = Z_1258655 x|_263 Z_288 = <x, y | x^1258655 = y^288 = 1, y x y^-1 = x^263>, of order 362,492,640 (1258655 = 5 x 19 x 13249), with connection set S = {x^286140 y^144, x^527543 y^2, x^579807 y^69, x^306882 y^14, x^402861 y^32, x^244866 y^149, x^463227 y^58}: every g in G is joined to g s and g s^-1 for s in S.

x^286140 y^144 is one of the 13,249 involutions in G, so the degree is 1 + 2 x 6 = 13. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 1258655 j is the group element x^i y^j (0 <= i < 1258655, 0 <= j < 288), so vertex 0 is the identity and the 288 runs of 1,258,655 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 288).

## How it was found

Small-complement group screen, bitmap simulated annealing with a movable involution in the one-involution Cayley family on Z_1258655 x|_263 Z_288. The search code is in [search/cayley_metacyclic_anneal](../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (2,356,202,160 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg13-dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search. It needs a few bits of memory per vertex, so `check` skips it unless asked:

    python3 tools/ddgraph.py check graphs/deg13-dia9 --max-vertices 362492640
