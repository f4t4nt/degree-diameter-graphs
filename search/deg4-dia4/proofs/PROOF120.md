# The reduced 120 puzzle has no solution: a short proof

This proves by hand what `design120_reduced.py` (../sat_closure/) checks by computer (section 11 of
`TEMPLATE_CALCULUS.md`): the four partition families (A)-(D) cannot be satisfied for any choice of
the OA bijections. The proof is local (it looks at one cell moved by kappa), uses all four families
but at most nine of the 512 exactly-one constraints, and uses nothing about kappa except that it is
not the identity. Every step was re-checked by `verify_proof.py` (../sat_closure/proof120/) for the three phi13
representatives and all 6912 phi23 (section 6).

## 1. Statement

Let P be the 16 cells of a 4x4 grid. A *star* is a row or a column; a *perfect matching* (PM) is a set
of four cells meeting every star exactly once. A bijection phi: P -> P' between two grids is an *OA
bijection* if it maps every star of P onto a PM of P'.

Fix OA bijections phi12: P1 -> P2, phi13: P1 -> P3, phi23: P2 -> P3 (all three port sets are copies
of the grid), put kappa = phi13^-1 phi23 phi12, a permutation of P1, and call the eight sets
M = phi12^-1(S), S a star of P2, the *OA matchings*.

**Theorem.** There are no permutations alpha, rho, sigma of P1 such that for every star S of P1 and
every OA matching M

    (A)  S, rho(S), alpha^-1(S), sigma(S)             partition P1,
    (B)  S, alpha(S), kappa rho(S), kappa sigma(S)     partition P1,
    (C)  M, alpha(M), rho^-1(M), sigma^-1(M)           partition P1,
    (D)  kappa(M), alpha(M), kappa rho(M), sigma^-1(M) partition P1.

More precisely, (A)-(D) force kappa to be the identity (Proposition, section 4), and kappa is never
the identity for OA bijections (Lemma 3).

## 2. Notation and the facts used

For a cell x of P1 write R(x), C(x) for its row and column, and L1(x), L2(x) for the two OA matchings
containing x: L1(x) = phi12^-1(row of phi12(x)), L2(x) = phi12^-1(column of phi12(x)). In Latin
square language L1(x), L2(x) are the symbol classes of x in the two orthogonal Latin squares behind
phi12. Call these four sets the *lines* of x.

**(F0)** phi is an OA bijection iff phi^-1 is: both say |phi(S) ∩ T| = 1 for every star S of P and
every star T of P'. In particular the OA matchings are PMs of P1.

**(F1)** Rows, columns, L1-matchings and L2-matchings are four partitions of P1 into 4-sets. So every
cell has exactly one line in each class, and y ∈ L_i(x) iff L_i(y) = L_i(x) (likewise for R, C).

**(F2)** Two lines of x from different classes meet only in x: R(x) ∩ C(x) = {x} trivially; a PM
meets a star in one cell, which handles L_i(x) against R(x) and C(x); and
L1(x) ∩ L2(x) = phi12^-1(row ∩ column of phi12(x)) = {x}. Hence a cell z != x lies on at most one
line of x.

### Reading the families one cell at a time

"X1, X2, X3, X4 partition P1" for four 4-sets says that every cell z lies in exactly one Xi, and
z ∈ pi(X) iff pi^-1(z) ∈ X. So (A)-(D) say, for every star S, OA matching M and cell z:

    (A)[S; z]  exactly one of   z ∈ S,            rho^-1(z) ∈ S,    alpha(z) ∈ S,               sigma^-1(z) ∈ S
    (B)[S; z]  exactly one of   z ∈ S,            alpha^-1(z) ∈ S,  rho^-1(kappa^-1 z) ∈ S,     sigma^-1(kappa^-1 z) ∈ S
    (C)[M; z]  exactly one of   z ∈ M,            alpha^-1(z) ∈ M,  rho(z) ∈ M,                 sigma(z) ∈ M
    (D)[M; z]  exactly one of   kappa^-1(z) ∈ M,  alpha^-1(z) ∈ M,  rho^-1(kappa^-1 z) ∈ M,     sigma(z) ∈ M

This is literally what `design120_reduced.solve_reduced` encodes. When z lies in S, respectively
in M, the first alternative holds and the other three must fail:

    (A0)  z ∈ S  =>  rho^-1(z), alpha(z), sigma^-1(z) are not in S,
    (C0)  z ∈ M  =>  alpha^-1(z), rho(z), sigma(z) are not in M.

## 3. Three lemmas

**Lemma 1 (alpha follows kappa; uses (A), (B)).** If z = kappa(x) lies in neither R(x) nor C(x),
then alpha(x) = z.

*Proof.* Apply (B)[R(x); z]. Its first alternative fails (z ∉ R(x)), and since kappa^-1(z) = x the
third and fourth read rho^-1(x) ∈ R(x) and sigma^-1(x) ∈ R(x), which fail by (A0) at x with
S = R(x). Hence alpha^-1(z) ∈ R(x). The same with C(x) gives alpha^-1(z) ∈ C(x), so
alpha^-1(z) ∈ R(x) ∩ C(x) = {x} by (F2). ∎

**Lemma 2 (rho^-1 follows kappa; uses (C), (D)).** If z = kappa(x) lies in neither L1(x) nor L2(x),
then rho^-1(x) = z.

*Proof.* By (F1), z ∉ L_i(x) means x ∉ L_i(z), i.e. kappa^-1(z) ∉ L_i(z). Apply (D)[L_i(z); z] for
i = 1, 2. Its first alternative fails, and the second and fourth, alpha^-1(z) ∈ L_i(z) and
sigma(z) ∈ L_i(z), fail by (C0) at z with M = L_i(z). Hence rho^-1(x) = rho^-1(kappa^-1 z) lies in
both L1(z) and L2(z), whose intersection is {z} by (F2). ∎

**Lemma 3 (a composite of two OA bijections is not an OA bijection).** phi13 phi12^-1: P2 -> P3 is
not an OA bijection. Equivalently, kappa != id for every OA bijection phi23.

*Proof.* kappa = id means phi23 = phi13 phi12^-1, so the two statements are equivalent. Suppose
psi = phi13 phi12^-1 were an OA bijection. For stars S of P2 and T of P3, psi(S) is a PM, so
1 = |psi(S) ∩ T| = |phi12^-1(S) ∩ phi13^-1(T)|: every OA matching of phi12 meets every OA matching
of phi13 in exactly one cell. Now consider the six partitions of P1 into 4-sets

    rows,  columns,  L1 = {phi12^-1(rows of P2)},  L2 = {phi12^-1(columns of P2)},
                     L1' = {phi13^-1(rows of P3)}, L2' = {phi13^-1(columns of P3)}.

Two blocks taken from two different partitions of this list always meet in exactly one cell: a row
and a column trivially; a star and an L- or L'-block because those blocks are PMs (F0); L1 against L2,
and L1' against L2', as in (F2); L against L' by the assumption. Fix a cell x. The six blocks through
x each contain three further cells, and these six 3-sets are pairwise disjoint (a common cell y != x
would put two cells in the intersection of two blocks), so P1 \ {x} would contain 18 cells. It
contains 15. ∎

Lemma 3 is the classical bound of n+1 = 5 mutually orthogonal parallel classes on n^2 = 16 points,
i.e. at most n-1 = 3 mutually orthogonal Latin squares of order 4: the six classes would be rows,
columns and four MOLS.

## 4. Proof of the theorem

**Proposition.** If alpha, rho, sigma satisfy (A)-(D) for some permutation kappa of P1, then kappa is
the identity.

*Proof.* Let x be a cell, z = kappa(x), and suppose z != x. By (F2), z lies on at most one of the
four lines of x. Three cases.

*(a) z lies on an OA matching M = L_i(x) of x.* Then z is in neither R(x) nor C(x), so alpha(x) = z
by Lemma 1, i.e. alpha^-1(z) = x ∈ M. But z ∈ M, so (C0) says alpha^-1(z) ∉ M. Contradiction.

*(b) z lies on a star S of x (S = R(x) or C(x)).* Then z is in neither L1(x) nor L2(x), so
rho^-1(x) = z ∈ S by Lemma 2. But x ∈ S, so (A0) says rho^-1(x) ∉ S. Contradiction.

*(c) z lies on none of the four lines of x.* Lemmas 1 and 2 both apply: alpha(x) = rho^-1(x) = z.
Look at (A)[R(z); x]. Its first alternative fails (x ∉ R(z), since z ∉ R(x)), while the second and
third, rho^-1(x) ∈ R(z) and alpha(x) ∈ R(z), both hold. Two alternatives out of "exactly one".
Contradiction.

So kappa(x) = x for every cell x. ∎

*Proof of the theorem.* A solution of (A)-(D) would make kappa the identity by the Proposition, but
kappa = phi13^-1 phi23 phi12 is never the identity for OA bijections by Lemma 3. ∎

The theorem holds for every triple of OA bijections, not only for the orbit representatives; no
symmetry reduction is used.

## 5. What the proof used, and what it did not

- **All four families are needed, but only a few instances of each.** For one cell x with
  z = kappa(x) != x the proof reads five constraints in cases (a) and (b) and nine in case (c):
  (A)[R(x); x], (A)[C(x); x], (B)[R(x); z], (B)[C(x); z] for Lemma 1; (C)[L1(z); z], (C)[L2(z); z],
  (D)[L1(z); z], (D)[L2(z); z] for Lemma 2; and (C)[M; z], (A)[S; x] or (A)[R(z); x] for the final
  clash. (A0) and (C0) use the "at most one" half of the constraints, Lemmas 1-2 the "at least one"
  half. No family can be dropped altogether: there are phi23 for which (A)+(B)+(C) is satisfiable
  (32, 26 and 32 of the 6912 phi23 for the three phi13 representatives) and phi23 for which
  (A)+(C)+(D) is satisfiable (8, 34 and 32 of 6912); for the first representative every pair of
  families is satisfiable for some phi23, and (A)+(C) for all of them.
- **Nothing about kappa except kappa != id.** Lemmas 1, 2 and the Proposition hold for an arbitrary
  permutation kappa of P1. The OA structure of phi23 enters only through Lemma 3. This is sharp:
  with kappa = id the four families are satisfiable (`proof120/verify_general.py` prints a
  solution), so "kappa = id" is exactly what (A)-(D) say about kappa. In the graph, kappa = id reads
  phi13 = phi23 phi12: the three OA bijections would have to form a commuting triangle, and Lemma 3
  says that the composite of two OA bijections always sends some hub star to a set that hits a hub
  of the far block twice.
- **Only the incidence facts (F0)-(F2).** The proof never uses the fifth parallel class of the affine
  plane behind phi12 (the four common transversals of L1 and L2, which complete rows, columns, L1, L2
  to an affine plane of order 4; `explore1.py` lists them), nor that alpha, rho, sigma preserve it
  (they do: (A0) and (C0) applied at every cell push alpha(x), rho(x), sigma(x) off all four lines
  of x, hence onto the fifth line through x), nor any property of the three phi13 representatives
  (their planes are: the plane of phi12 with L1' = L1, L2' = L2; the other plane through rows and
  columns; the plane of phi12 with L1' = L1 and L2' the fifth class). Orthogonality of the two Latin
  squares behind phi12 is used only as L1(x) ∩ L2(x) = {x}, and inside Lemma 3.
- **Hub-freeness** (alpha, rho, sigma move every cell off its row and column) is (A0); (C0) is its
  twin for the OA matchings. Both are used only through the lines of the one cell x.

## 6. Computational checks

All scripts are in ../sat_closure/proof120/ and run from that directory with a Python that has
python-sat (CaDiCaL) installed; each prints its findings.

- `verify_proof.py {0,1,2}`: for each phi13 representative and all 6912 phi23,
  Lemma 3 (phi13 phi12^-1 is not among the 6912 OA bijections, so kappa != id), and for every cell x
  moved by kappa the exact constraint set of the applicable case (a), (b) or (c) is unsatisfiable by
  itself (103,680 moved cells per representative: 41,472 in case (a), 41,472 in (b), 20,736 in (c)).
  On the first 300 phi23 also Lemma 1 (the four (A)/(B) constraints force alpha(x) = kappa(x)) and
  Lemma 2 (the four (C)/(D) constraints force rho^-1(x) = kappa(x)) are checked as stated, by adding
  the negated conclusion.
- `verify_general.py`: (A)-(D) with kappa = id is satisfiable (solution printed);
  with kappa an arbitrary non-identity permutation it is unsatisfiable for all 120 transpositions,
  500 random permutations and 53 random 3-cycles; phi13 phi12^-1 is not an OA bijection for the
  fixed phi12 and all 6912 phi13 (by the Aut x Aut symmetry this is Lemma 3 for every pair); and
  the satisfiable counts for pairs of families (representative 0): (A)+(B) 776, (A)+(D) 200,
  (B)+(C) 200, (B)+(D) 200, (C)+(D) 48 of 6912, (A)+(C) all 6912.
- `triples.py ABC`, `triples.py ACD`: the satisfiable counts for the two
  family triples quoted in section 5.
- The exploration that led here (the affine-plane classes behind phi12, how the three phi13
  representatives sit relative to that plane, kappa types, and deletion-based minimal
  unsatisfiable subsets of 5 to 13 constraints) is not included; the 5-constraint cores are
  exactly case (a), and reading them off suggested Lemmas 1-2 and the case split.
