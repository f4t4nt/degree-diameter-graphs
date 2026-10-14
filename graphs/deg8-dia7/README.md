# Degree 8, Diameter 7

A graph with 138,050 vertices, maximum degree 8 and diameter 7. The degree-diameter table listed 137,247 vertices for this cell on October 14, 2026. The Moore bound for this cell is 1,098,057.

## Construction

The Cayley graph of G = Z_2761 x|_63 Z_50 = <x, y | x^2761 = y^50 = 1, y x y^-1 = x^63>, of order 138,050 (2761 = 11 x 251), with connection set S = {x^290 y^14, x^971 y^43, x^457 y^6, x^1728 y^16}: every g in G is joined to g s and g s^-1 for s in S.

S contains no involution (G has 11 involutions), so the degree is 2 x 4 = 8. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 2761 j is the group element x^i y^j (0 <= i < 2761, 0 <= j < 50), so vertex 0 is the identity and the 50 runs of 2,761 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 50).

## How it was found

Even-M metacyclic group with a tiny involution class chosen by a random-set screen, simulated annealing followed by exact single-slot completion; Cay(Z_2761 x|_63 Z_50) with four inverse pairs.

## Files

- `deg8-dia7-n138050.edges.txt.gz`: the 552,200 edges (gzipped), one "u v" per line in the numbering above, grouped by generator under "#" headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The edge list is checked to be exactly the construction in `construction.json`. The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

    python3 tools/ddgraph.py check graphs/deg8-dia7
