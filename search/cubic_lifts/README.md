# Cubic Lifts

I sweep regular voltage lifts of small cubic base graphs (multigraphs with loops, parallel edges and semi-edges on k
vertices) over every group of order n/k given by an explicit multiplication table, one representative per orbit of the
group's automorphisms, for cubic graphs of a target diameter; the sweep of the 5-vertex bases over the groups of order
40 gave the (3,7) graph on 200 vertices (a lift over Z_5 x| Z_8).

Files:

- `lift.c`: the evaluator and exhaustive enumerator (iterated stabilisers of Aut(G), lexicographically minimal
  representatives, loop voltages identified with their inverses) with `exhaust`, `sample` and `eval` modes; a lift vertex
  (u, g) has index u |G| + g.
- `groups.py`: every group of a given order as a multiplication table (cyclic, abelian, direct and semidirect products,
  dicyclic, small permutation and matrix groups), deduplicated by a canonical form of the Cayley-table graph, with Aut(G)
  as permutations, written to `groups/order_<m>/g<i>.txt` plus an index.
- `bases.py`: all connected cubic base pregraphs on k vertices up to isomorphism in lift's base format
  (`k nslots ntree`, then `u v type` lines, type 0 edge, 1 loop, 2 semi-edge; the first ntree slots are a spanning tree).
- `cubic_bases.py`: simple cubic base graphs on K vertices prefiltered by the non-backtracking walk condition (for the
  larger k of the diameter-6 and diameter-8 cells).
- `sweep.py`: the driver that runs `lift` over every (group, base) pair of the requested orders and appends JSON lines.

To run it yourself (seconds; the group and base enumerations are deterministic, so `groups/order_40/g12.txt` is
Z_5 x| Z_8 with |Z| = 4 and `bases/b5_0026.txt` is the base with edges 0-1, 0-2, 0-3, 1-4, 1-2, a semi-edge at 2 and
loops at 3 and 4):

```
pip install -r requirements.txt
gcc -O2 -o lift lift.c
python3 groups.py groups 40          # 14 groups of order 40 (the known number), groups/order_40/index.txt
python3 bases.py bases 5             # 43 connected cubic base pregraphs on 5 vertices
./lift groups/order_40/g12.txt bases/b5_0026.txt 7 exhaust 10 sol
# expected: RESULT mode=exhaust m=40 k=5 n=200 D=7 naut=80 ninv=1 nfree=4 evaluated=487 solutions=2 best_defect=0 ...
#           BEST defect=0 voltages: 1 4 10 18          (sol_1.edges and sol_2.edges: 200 vertices, cubic, diameter 7)
```

Group elements are numbered by `groups.py` (element 0 is the identity; the (3,7) record uses the index a^i x^h -> 8 i + h
of Z_5 x| Z_8). A whole-order sweep is `python3 sweep.py --D 7 --n 200 --kmax 8 --groups groups --bases bases`, which
needs the groups of every order dividing 200 and the bases for every k.
