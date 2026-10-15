# Degree 16, Diameter 9

A graph with 2,343,053,633 vertices, maximum degree 16 and diameter 9. The degree-diameter table listed 2,120,003,902 vertices for this cell on October 14, 2026. The Moore bound for this cell is 43,935,267,857.

## Construction

The Cayley graph of G = Z_5060591 x|_949 Z_463 = <x, y | x^5060591 = y^463 = 1, y x y^-1 = x^949>, of order 2,343,053,633 (5060591 is prime), with connection set S = {x^677132 y^226, x^3285509 y^90, x^3540687 y^80, x^40693 y^78, x^1582586 y^22, x^3137470 y^103, x^3113179 y^174, x^4144856 y^377}: every g in G is joined to g s and g s^-1 for s in S.

G has odd order, so it has no involutions, and the degree is 2 x 8 = 16. The graph is vertex-transitive: G acts on it by left multiplication.

Vertex i + 5060591 j is the group element x^i y^j (0 <= i < 5060591, 0 <= j < 463), so vertex 0 is the identity and the 463 runs of 5,060,591 consecutive numbers are the cosets of <x>; a generator x^a y^b joins run j to run j + b (mod 463).

## How it was found

Odd-order group with small complement M, bitmap simulated annealing with targeted moves on Z_5060591 x|_949 Z_463. The search code is in [search/cayley_metacyclic_anneal](../../search/cayley_metacyclic_anneal).

## Files

- No edge list is stored (18,744,429,064 edges); the construction defines it, and `python3 tools/ddgraph.py edges graphs/deg16-dia9` streams it in the numbering above.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

The graph is vertex-transitive, so the eccentricity of the identity is the diameter: a breadth-first search over the graph defined by `construction.json` decides it.

`verification.txt` is the output of that search. It needs a few bits of memory per vertex, so `check` skips it unless asked:

    python3 tools/ddgraph.py check graphs/deg16-dia9 --max-vertices 2343053633
