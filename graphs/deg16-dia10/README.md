# Degree 16, Diameter 10

A graph with 29,899,503,497 vertices, maximum degree 16 and diameter 10. The degree-diameter table listed 22,500,000,906 vertices for this cell on October 14, 2026. The Moore bound for this cell is 659,029,017,857.

## Construction

The Cayley graph of G = Z_24170981 x|_16573517 Z_1237 = <x, y | x^24170981 = y^1237 = 1, y x y^-1 = x^16573517>, of order 29,899,503,497 (24170981 is prime), with connection set S = {x^10866025 y^532, x^21248118 y^539, x^13248079 y^918, x^21841737 y^1181, x^1620224 y^1059, x^12379694 y^1095, x^17981217 y^470, x^3158481 y^1021}: every g in G is joined to g s and g s^-1 for s in S.

G has odd order, so it has no involutions, and the degree is 2 x 8 = 16. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 24170981 j is the group element x^i y^j (0 <= i < 24170981, 0 <= j < 1237), so vertex 0 is the identity and the 1237 runs of 24,170,981 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 1237).

## How it was found

Odd-order group with prime N and M and M about 0.40 n^(1/3), two-plane bitmap simulated annealing with targeted moves on Z_24170981 x|_16573517 Z_1237. The search code is in [search/cayley_metacyclic_anneal](../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (239,196,027,976 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg16-dia10` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search. It needs a few bits of memory per vertex, so `check` skips it unless asked:

    python3 tools/ddgraph.py check graphs/deg16-dia10 --max-vertices 29899503497
