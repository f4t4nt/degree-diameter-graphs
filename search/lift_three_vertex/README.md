# Lift Three Vertex

I extend the two-vertex lift sweeps to 3-vertex base graphs: every connected cubic-like base shape on three vertices
(links with multiplicities, loops and semi-edges) is enumerated up to isomorphism, turned into a `bitlift` base spec, and
its voltage assignments over a metacyclic group G(N, M, r, t) are swept with the same sequential-stabiliser reduction and
the same evaluator (`../lift_two_vertex/bitlift.c`, compiled with `-DMAXK=4`). The degree-4 base with a triple link, a
path and a loop over the one-involution group Z_73 x|_18 Z_36 gave the (4,9) record on 7884 vertices.

Files:

- `shapes.py`: the 61 connected 3-vertex degree-5 base pregraphs up to isomorphism with their automorphism group order,
  free-slot structure and base spec (`--json` for machine use); degree-4 shapes are written by hand as specs.
- `gen_k3.py`: reduced candidate prefixes for any base spec (group automorphisms, loop inversion, permutations of
  parallel slots; shape automorphisms are not used, so the family is over-covered by at most |Aut(shape)|).
- `check_k3.py`: brute-force validation of the reduction on a small group (group axioms, every symmetry a lift
  isomorphism, orbit coverage).
- `make_edges_k3.py`: independent edge-list builder for a k-vertex lift (lift vertex (u, g) has index u |G| + index(g)).

Base spec of the (4,9) record: `k=3;deg=4;loop:2;link*:0,1;link*:0,2;link:1,2;link:0,1;link:0,1` (free slots in order:
the loop at vertex 2, the link 1-2, the two further links 0-1; the starred links are the spanning tree with trivial
voltage). The full family over Z_73 x|_18 Z_36 has about 1.8e9 candidates after reduction; the control below sweeps the
sub-block with the loop fixed to y^3 (index 219) and the link 1-2 fixed to x y^6 (index 439).

To run it yourself (about half a minute with AVX-512; re-finds the record assignment):

```
pip install -r requirements.txt
gcc -O3 -march=native -DNW=2 -DMAXK=4 -o bitlift_nw2 ../lift_two_vertex/bitlift.c -lm
python3 -c "
for c in range(2628): print(219, 439, c, 0, 2628)" > pref_b3g.txt
./bitlift_nw2 --meta 73 36 18 --base "k=3;deg=4;loop:2;link*:0,1;link*:0,2;link:1,2;link:0,1;link:0,1" --D 9 --mode prefix --prefix pref_b3g.txt
# expected: FOUND diam<=9 (defect 0) group=Z73 x|_18 Z36 ... voltages=[219,0,0,439,779,1062]
#           FOUND ... voltages=[219,0,0,439,1062,779]        (the same lift with the two parallel links swapped)
#           EXHAUSTIVE ... candidates=6903756 ... best_defect=0 ... found=2 complete=1
python3 make_edges_k3.py 73 36 18 0 "k=3;deg=4;loop:2;link*:0,1;link*:0,2;link:1,2;link:0,1;link:0,1" "[219,0,0,439,779,1062]" 7884.edges
```

The reduced prefixes of a whole family come from `python3 gen_k3.py N M r t "SPEC" --out pref.txt` (add `--count-only`
to size it first); `python3 shapes.py` lists the degree-5 shapes.
