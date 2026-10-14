# Degree 5, Diameter 7

A graph with 6,000 vertices, maximum degree 5 and diameter 7. The degree-diameter table listed 5,896 vertices for this cell on October 14, 2026. The Moore bound for this cell is 27,306.

## Construction

A regular voltage lift over the group G described below, of order 3,000. The base graph has vertices 0 and 1, with three parallel links 0 -> 1 (voltages 1, v(0,1) c z y^4, v(1,4) c z^4 y^7), a loop at 0 (voltage y) and a loop at 1 (voltage v(0,2) c z^4 y^4).

The lift has vertex set {0, 1} x G and joins (u, g) to (v, g a) for every dart u -> v with voltage a; a loop with voltage a joins (v, g) to (v, g a) and (v, g a^-1). G acts on each fibre {v} x G by left multiplication.

The group is G = ((F_25 x| <c>) x <z>) x| <y> of order 3000: F_25 = F_5^2 as column vectors, c of order 3 acting on F_25 by B = [[0,4],[1,4]] (so F_25 x| <c> is a Frobenius group of order 75), z central of order 5 in K = (F_25 x| <c>) x <z>, and y of order 8 acting on K by (w, c^i, z^j) -> (A w, c^(2i), z^(2j)) with A = [[1,1],[2,4]]. G is not metacyclic; it has trivial centre and 25 involutions. The word v(a,b) c^i z^j y^m stands for the element ((a,b), c^i, z^j) y^m.

Vertex 3000 v + index(g) is the pair (v, g), where g = ((a,b), c^i, z^j) y^m has index (((5 a + b) 3 + i) 5 + j) 8 + m, so vertex 3000 v is (v, identity).

## How it was found

Complete sweep of all shape-103 two-vertex lift assignments (one loop and three links per base vertex, one link trivial) on the non-metacyclic group ((F_25 x| C_3) x Z_5) x| C_8 up to symmetry, the group selected by the conjugacy-class-count (commuting probability) rule.

## Files

- `deg5-dia7-n6000.edges.txt`: the 15,000 edges, one "u v" per line in the numbering above, grouped by voltage under "#" headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter; the edge list is also checked to be exactly the construction in `construction.json`.

    python3 tools/ddgraph.py check graphs/deg5-dia7
