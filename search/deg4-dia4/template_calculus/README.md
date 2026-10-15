# Template Calculus (Counting Side)

I derive what a degree-4 diameter-4 graph built from three subdivided $`K_{4,4}`$ blocks must look like,
from a reach identity and a port lemma down to 54 skeletons, and these scripts carry out the
finite checks of `../proofs/TEMPLATE_CALCULUS.md` sections 1-8, 11b and the two-block analysis.

- `verify.py`: independent all-pairs BFS verifier plus the three-block decomposition of a graph,
  its connector profiles and the hub-to-port reach identity (`verify.py corpus` applies it to the 36
  classes; `verify.py edges FILE...` to any edge list). Other scripts import it.
- `counting_ilp.py`: the counting relaxation over connector profiles (integer program with the
  per-block slot sums, the budget rows and hub-to-hub coverage); gives budget max 768 at 105.
- `shapes105.py`: the 84 special-connector configurations at 105 that survive the exact budget rows
  and per-hub compensation; `port_lemma.py` reduces them to the 21 triangle configurations.
- `surgery.py`: pair-to-triangle surgery on the 104s (`t1`: every pair of every class);
  `unit_structure.py`: the twin structure predicted by the structure theorem, checked on all 36.
- `unit_designs2.py`: orderly enumeration of tight designs (three fixed-point-free involutions
  on 16 units with the side and label conditions) under the true label group, nauty-reduced: 54
  skeletons (`../sat_closure/data/skeletons54.json`); `known_designs2.py` maps the 36 classes to
  their skeletons (50: 8 classes, 51: 18, 53: 10).
- `blocks_relax.py`: the free blocks (subdivided $`K_5`$, octahedron, $`K_{4,4}`$) and the hub-reach
  relaxation over block multisets. `potential3.py table`, `check_potential_table.py`: the potential
  table of section 13.2 (every connector-edge potential is $`\le 0`$, with equality exactly for
  complementary triple pairs and cross-double to same-block-double edges), with exact rationals.
- `twoblock/`: two $`K_{4,4}`$ blocks plus an arbitrary remainder: the two-block reach identity and
  relaxation (`twoblock_ilp.py`, feasible up to $`n = 142`$, no budget cut), calibration on the 36
  classes (`known104_twoblock.py`), a free-remainder SAT model with controls (`twoblock_sat.py`,
  `fixedH_control.py`, `negative_control.py`, `plan_c_105.py`, `analyze_solution.py`); their
  docstrings cite sections of a separate two-block note that is not staged; its conclusion is the
  one above.

Work from this directory; the scripts print JSON or write to `runs/` (`counting_ilp.py`,
`shapes105.py` and `blocks_relax.py` need scipy, the rest pynauty, networkx and numpy). To run it yourself:

```
python verify.py corpus > verify36.json            # 36 graphs: diameter 4, identity 768 = 768
python counting_ilp.py                             # relaxation: max budget 768 at n = 105
mkdir -p runs && python shapes105.py > runs/shapes105_v3.json && python port_lemma.py   # 84 -> 21
python unit_designs2.py > skeletons.json           # 54 skeletons
python known_designs2.py                           # classes -> skeletons 50, 51, 53
```
