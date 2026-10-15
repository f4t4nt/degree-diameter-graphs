# Degree 4, Diameter 4

A graph with 104 vertices, maximum degree 4 and diameter 4. This is my entry in the degree-diameter table for this cell (104 vertices on October 14, 2026). The Moore bound for this cell is 161.

## Construction

Three disjoint copies (blocks A, B and C) of the subdivided complete bipartite graph K_4,4, and 32 connector vertices in 16 matched pairs. In a block, each of the hubs L0..L3 and R0..R3 meets four ports, and port P_xy subdivides the edge L_x R_y. Every port also meets two connectors, and every connector meets its partner and three ports.

Block A is 0..23 (L0..L3 = 0..3, R0..R3 = 4..7, port P_xy = 8 + 4x + y), block B is 24..47 and block C is 48..71 in the same pattern, and connector pair k (0..15) is 72 + 2k and 73 + 2k.

This graph is class 29 of the 36 graphs of degree 4 and diameter 4 on 104 vertices built from three subdivided K_4,4 blocks with no port-port edge, one of the two classes with an automorphism group of order 32. All 36 are in [search/deg4-dia4/classes](../../search/deg4-dia4/classes), and [search/deg4-dia4](../../search/deg4-dia4) has the proofs that this architecture admits no graph on 105, 106, 107 or 108 vertices.

## How it was found

Symmetric Z_8 voltage lift over a 13-vertex base pregraph with defect-targeted moves and basin hopping, run until a seed stalls at a small defect, then symmetry-breaking 2-opt edge repair closes the last pairs; the 36 isomorphism classes come from enumerating every diameter-4 realisation (SAT) of the three tight three-block skeletons that admit one and sorting them by canonical form. The search code is in [search/deg4-dia4](../../search/deg4-dia4).

## Files

- `deg4-dia4-n104.edges.txt`: the 208 edges, one "u v" per line in the numbering above, grouped by edge type under "#" headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter.

    python3 tools/ddgraph.py check graphs/deg4-dia4
