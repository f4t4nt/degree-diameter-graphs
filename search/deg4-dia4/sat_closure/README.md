# SAT Closure (Completion Encodings, Shape Enumerations, Impossibility Checks)

I decide each finite case that the template calculus leaves open with a SAT encoding whose every
clause is a walk-of-length-at-most-4 condition, validated on the 36 classes at 104 (SAT with their
own wiring, SAT breaking after one port swap) and on negative controls; the statements proved are in
`../proofs/TEMPLATE_CALCULUS.md`, sections 9 to 13. Work from this directory (outputs go to
`runs/`); the scripts need `python-sat` with CaDiCaL and `pynauty` (`../requirements.txt`).

Files:

- `sat105.py MODE [SKELETONS]` (modes `104`, `105`, `120`; `sat106.py` adds `106`..`119`): completion
  of a skeleton (`data/skeletons54.json`) over the 10 splits per unit and the 15 triangle-gadget
  options. `104 all` is SAT exactly on skeletons 50, 51, 53; `105`, `106`, `107`, `108` and `120 all`
  are UNSAT on all 54. `validate.py`: the encoder on the 36 classes plus a gadget control.
- `enumerate104.py K`, `classify104.py K`: every diameter-4 realisation of skeleton $`K`$ reduced to
  isomorphism classes with nauty: 8 + 18 + 10 = 36 classes, the files of `../classes/`.
- `typed_slack_sat.py`: the typed encoder for shapes with budget slack (gadget profiles and edges
  fixed, unit types free, redundant port-lemma and exactness clauses, tight mode at budget 768), with
  `validate` and `control` modes. `typed_sat_shapes.py SHAPES --out=runs/X.json --cap=S` applies it
  to a shape list: every shape of `data/shapes106.json` (54), `shapes107.json` (729),
  `shapes108_776plus.json` (3,328), `shapes108_772.json.gz` (2,325) and `shapes108_768.json.gz`
  (11,163) is UNSAT, which together with the skeleton completions closes 105 to 108.
  `cube_sat108e.py` is the case-split fallback for the hardest shape; `merge108e.py` merges the
  verdict files of the batches per shape.
- Shape enumerations behind those lists: `components107.py` (typed connector components by
  potential), `shapes106.py`, `shapes107.py`, `shapes_fast.py`, `shapes108.py` (units plus gadget
  multisets under the budget rows and the port-lemma matching; `recheck_tight.py` drops the unsound
  assumption (c)); for 108 at budgets 772 and 768 the window pipeline `windows108.py`,
  `export_windows.py` -> `clusters108.c` (7,702,003 clusters of cost $`\le 12`$) -> `combine108.c`,
  `combine_groups.c`, `combine_classes.c` -> `c2jsonl.py`, `shapes_from_c.py` -> `windows108.py shapes`,
  each stage validated against the previous one and against the 106 list.
- `design120.py`, `design120_reduced.py`: the 120-vertex all-cross-double design, UNSAT for the three
  $`\varphi_{13}`$ representatives and in all 20,736 reduced instances (`../proofs/PROOF120.md`; checks in
  `proof120/`). `oct4/`: four octahedra at 108 (potential table, ILP cross-check, holder table, typed
  SAT over the 11 skeletons, all UNSAT). `design108.py`, `design110.py`: a profile-free exactness
  model of the four-octahedron design at 108 and a structured model of an octahedron plus two
  $`K_{4,4}`$ blocks at 110, neither part of a proof chain; `direct_sat.py`: the unstructured
  encoding, included to show why every decision is structured. `check_coholding.py`: none of the 36
  classes has a connector pair co-holding two ports.

To run it yourself:

```
python sat105.py 104 50,51,53 && python sat105.py 105 all && python validate.py
# expected: SAT x3, UNSAT x54, then the 36 controls
python typed_sat_shapes.py data/shapes106.json --out=runs/typed106.json --cap=900
# expected: 54 x UNSAT
python design120.py && python design120_reduced.py && (cd oct4 && python typed_sat.py oct4 all)
# expected: UNSAT throughout
```
