# Degree 10, Diameter 10

A graph with 265,612,889 vertices, maximum degree 10 and diameter 10. The degree-diameter table listed 239,005,058 vertices for this cell on October 14, 2026. The Moore bound for this cell is 4,358,480,501.

## Construction

The Cayley graph of G = Z_1111351 x|_8147 Z_239 = <x, y | x^1111351 = y^239 = 1, y x y^-1 = x^8147>, of order 265,612,889 (1111351 is prime), with connection set S = {x^648877 y^104, x^1110715 y^82, x^841493 y^123, x^484351 y^59, x^482973 y^144}: every g in G is joined to g s and g s^-1 for s in S.

G has odd order, so it has no involutions, and the degree is 2 x 5 = 10. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 1111351 j is the group element x^i y^j (0 <= i < 1111351, 0 <= j < 239), so vertex 0 is the identity and the 239 runs of 1,111,351 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 239).

## How it was found

Odd-order group with small complement M, bitmap simulated annealing with targeted moves on Z_1111351 x|_8147 Z_239. The search code is in [search/cayley_metacyclic_anneal](../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (1,328,064,445 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg10-dia10` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search. It needs a few bits of memory per vertex, so `check` skips it unless asked:

    python3 tools/ddgraph.py check graphs/deg10-dia10 --max-vertices 265612889
