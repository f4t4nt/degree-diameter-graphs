# Proofs

I wrote up the theory behind the 104 classification and the 105 to 108 closures in three notes;
the scripts they cite are in `../template_calculus/` and `../sat_closure/`.

- `TEMPLATE_CALCULUS.md`: the template calculus for the three-K4,4-block family. Sections 1-10
  derive the hub-to-port identity, the counting relaxation, the port lemma, the structure theorem
  (54 skeletons) and the completion SAT, and prove that no 105-vertex member exists. Section 11
  settles the 120-vertex all-cross-double shape, 11b the four-octahedron design at 108, and section
  13 collects the closures at 106, 107 and 108 (potential function, component and window
  enumerations, typed SAT, parity of the budget).
- `PROOF120.md`: the hand proof that the reduced 120 puzzle has no solution: the four partition
  families force the three orthogonal-array bijections to commute, and a composite of two OA
  bijections is never an OA bijection because that would need four MOLS of order 4.
- `SOUNDNESS_A_PRIME.md`: why the refined excess matching used to filter the 108 shape lists is a
  necessary condition for every member at every budget.

The notes are kept as written during the work; only machine-specific details were removed.
