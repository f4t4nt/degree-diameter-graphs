# Voltage Lift Plus 2-Opt Repair (The Search That Found 104)

I found the first 104-vertex graphs of degree 4 and diameter 4 by a symmetric voltage-lift search
that stalls at a small defect, followed by a symmetry-breaking 2-opt repair that closes the last
pairs. Neither half works alone: the symmetric search never reaches defect 0 at 104, and plain
local search from scratch walls far below 100.

- `lift.c`: voltage-lift search. A base pregraph on B vertices with link, loop and semi-edge
  classes carries voltages in a finite group (cyclic, dihedral, dicyclic, Frobenius, products);
  the lift on B x |G| vertices is 4-regular when every base vertex has degree 4 counting a loop
  twice. The defect (ordered pairs at distance > 4) is evaluated from orbit representatives with
  bitset ball expansion; moves are defect-targeted voltage changes, an exhaustive polish over all
  class voltages and endpoint swaps, link/semi-edge conversions, and basin hopping from the best
  state. Usage: `./lift GROUP B F SEED SECONDS CKPT` (`cyc:8 13 0` is the record geometry: Z_8
  over 13 base vertices, no fixed vertices); `SECONDS = 0` dumps the checkpoint's graph.
- `repair.c`: plain 2-opt / 3-opt descent on the full graph with no symmetry, targeted at the
  vertices that still have far pairs, with kicks from the best state. Usage:
  `./repair EDGES SEED SECONDS BEST` (extra modes: `a` anneal, `g` exhaustive single-vertex
  insertion, `2` double insertion, `w` plateau walk).
- `z8_b13_defect16.ck`: the recorded stall, a Z_8 lift of a 13-vertex base with 28 classes at
  defect 16 (8 unordered pairs, all inside one fibre); the repair closes it in a fraction of a
  second from every seed tried, landing in classes 0, 9, 10 or 11 of `../classes/`.
- `run_control.sh`: compiles everything, dumps the stall, repairs it and verifies the result with
  the repository's all-pairs BFS verifier, `tools/verify/verify_bfs.c`.

The lift search is a stochastic basin-hopping search; a fresh run of `./lift cyc:8 13 0 SEED 60 x.ck`
reaches defects of 30 to 180 for most seeds within a minute and stalls at defect <= 24, where the
repair can take over, for roughly one seed in eighty. The seed that produced the recorded stall does
not reproduce it with this source (the search was edited after that run), which is why the stall
itself is included.

To run it yourself:

```
sh run_control.sh            # dump the stall (defect 16), repair to 0, verify by BFS
sh run_control.sh 5          # a different repair seed; expect "unordered pairs at distance > 4: 0"
./lift cyc:8 13 0 123 60 s123.ck    # a fresh one-minute lift search from seed 123
./repair s123.ck.edges 1 60 out.edges   # worth trying when the lift stalled at defect <= 24
```
