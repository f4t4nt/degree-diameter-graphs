# The reduced 120 puzzle has no solution: a short proof

This proves by hand what `design120_reduced.py` (../sat_closure/) checks by computer (section 11 of
`TEMPLATE_CALCULUS.md`): the four partition families (A)-(D) cannot be satisfied for any choice of
the OA bijections. The proof is local (it looks at one cell moved by $`\kappa`$), uses all four families
but at most nine of the 512 exactly-one constraints, and uses nothing about $`\kappa`$ except that it is
not the identity. Every step is re-checked by `verify_proof.py` (../sat_closure/proof120/) for the
three $`\varphi_{13}`$ representatives and all 6912 $`\varphi_{23}`$ (section 6).

## 1. Statement

Let $`P`$ be the set of 16 cells of a $`4 \times 4`$ grid. A *star* is a row or a column of the grid; a *perfect
matching* (PM) is a set of four cells meeting every star exactly once. A bijection $`\varphi: P \to P'`$ between two
such grids is an *OA bijection* if it maps every star of $`P`$ onto a PM of $`P'`$.

Fix OA bijections $`\varphi_{12}: P_1 \to P_2`$, $`\varphi_{13}: P_1 \to P_3`$, $`\varphi_{23}: P_2 \to P_3`$ (all
three port sets are copies of the grid), put $`\kappa = \varphi_{13}^{-1} \varphi_{23} \varphi_{12}`$, a permutation
of $`P_1`$, and call the eight sets $`M = \varphi_{12}^{-1}(S)`$, $`S`$ a star of $`P_2`$, the *OA matchings*.

**Theorem.** There are no permutations $`\alpha, \rho, \sigma`$ of $`P_1`$ such that for every star $`S`$ of
$`P_1`$ and every OA matching $`M`$

```math
\begin{array}{ll}
\text{(A)} & S,\ \rho(S),\ \alpha^{-1}(S),\ \sigma(S) \ \text{ partition } P_1, \\
\text{(B)} & S,\ \alpha(S),\ \kappa\rho(S),\ \kappa\sigma(S) \ \text{ partition } P_1, \\
\text{(C)} & M,\ \alpha(M),\ \rho^{-1}(M),\ \sigma^{-1}(M) \ \text{ partition } P_1, \\
\text{(D)} & \kappa(M),\ \alpha(M),\ \kappa\rho(M),\ \sigma^{-1}(M) \ \text{ partition } P_1.
\end{array}
```

More precisely, the families (A)-(D) force $`\kappa`$ to be the identity (Proposition, section 4), and $`\kappa`$
is never the identity for OA bijections (Lemma 3, section 3).

## 2. Notation and the facts used

For a cell $`x`$ of $`P_1`$ write $`R(x), C(x)`$ for its row and column, and $`L_1(x), L_2(x)`$ for the two OA
matchings containing $`x`$: $`L_1(x) = \varphi_{12}^{-1}(\text{row of } \varphi_{12}(x))`$,
$`L_2(x) = \varphi_{12}^{-1}(\text{column of } \varphi_{12}(x))`$. In Latin square language $`L_1(x), L_2(x)`$ are
the symbol classes of $`x`$ in the two orthogonal Latin squares behind $`\varphi_{12}`$. Call these four sets
the *lines* of $`x`$.

**(F0)** $`\varphi`$ is an OA bijection iff $`\varphi^{-1}`$ is: both say $`|\varphi(S) \cap T| = 1`$ for every star
$`S`$ of $`P`$ and every star $`T`$ of $`P'`$. In particular the OA matchings are PMs of $`P_1`$.

**(F1)** Rows, columns, $`L_1`$-matchings and $`L_2`$-matchings are four partitions of $`P_1`$ into 4-sets. So
every cell has exactly one line in each class, and $`y \in L_i(x)`$ iff $`L_i(y) = L_i(x)`$ (likewise for
$`R, C`$).

**(F2)** Two lines of $`x`$ from different classes meet only in $`x`$: $`R(x) \cap C(x) = \{x\}`$ trivially; a PM
meets a star in one cell, which handles $`L_i(x)`$ against $`R(x)`$ and $`C(x)`$; and
$`L_1(x) \cap L_2(x) = \varphi_{12}^{-1}(\text{row} \cap \text{column of } \varphi_{12}(x)) = \{x\}`$. Hence a cell
$`z \ne x`$ lies on at most one line of $`x`$.

### Reading the families one cell at a time

"$`X_1, X_2, X_3, X_4`$ partition $`P_1`$" for four 4-sets says that every cell $`z`$ lies in exactly one $`X_i`$,
and $`z \in \pi(X)`$ iff $`\pi^{-1}(z) \in X`$. So (A)-(D) say, for every star $`S`$, OA matching $`M`$ and cell $`z`$:

```math
\begin{array}{llllll}
\text{(A)}[S; z] & \text{exactly one of} & z \in S, & \rho^{-1}(z) \in S, & \alpha(z) \in S, & \sigma^{-1}(z) \in S \\
\text{(B)}[S; z] & \text{exactly one of} & z \in S, & \alpha^{-1}(z) \in S, & \rho^{-1}(\kappa^{-1} z) \in S, & \sigma^{-1}(\kappa^{-1} z) \in S \\
\text{(C)}[M; z] & \text{exactly one of} & z \in M, & \alpha^{-1}(z) \in M, & \rho(z) \in M, & \sigma(z) \in M \\
\text{(D)}[M; z] & \text{exactly one of} & \kappa^{-1}(z) \in M, & \alpha^{-1}(z) \in M, & \rho^{-1}(\kappa^{-1} z) \in M, & \sigma(z) \in M
\end{array}
```

This is literally what `design120_reduced.solve_reduced` encodes. When $`z`$ lies in $`S`$, respectively
in $`M`$, the alternative $`z \in S`$ (respectively $`z \in M`$) holds and the other three must fail:

```math
\begin{array}{ll}
\text{(A0)} & z \in S \ \Rightarrow\ \rho^{-1}(z),\ \alpha(z),\ \sigma^{-1}(z) \text{ are not in } S, \\
\text{(C0)} & z \in M \ \Rightarrow\ \alpha^{-1}(z),\ \rho(z),\ \sigma(z) \text{ are not in } M.
\end{array}
```

## 3. Three lemmas

**Lemma 1 ($`\alpha`$ follows $`\kappa`$; uses (A), (B)).** If $`z = \kappa(x)`$ lies in neither $`R(x)`$ nor $`C(x)`$,
then $`\alpha(x) = z`$.

*Proof.* Apply (B)$`[R(x); z]`$. Its alternative $`z \in R(x)`$ fails, and since $`\kappa^{-1}(z) = x`$
the third and fourth read $`\rho^{-1}(x) \in R(x)`$ and $`\sigma^{-1}(x) \in R(x)`$, which fail by (A0) at $`x`$ with
$`S = R(x)`$. Hence $`\alpha^{-1}(z) \in R(x)`$. The same with $`C(x)`$ gives $`\alpha^{-1}(z) \in C(x)`$, so
$`\alpha^{-1}(z) \in R(x) \cap C(x) = \{x\}`$ by (F2). $`\blacksquare`$

**Lemma 2 ($`\rho^{-1}`$ follows $`\kappa`$; uses (C), (D)).** If $`z = \kappa(x)`$ lies in neither $`L_1(x)`$ nor
$`L_2(x)`$, then $`\rho^{-1}(x) = z`$.

*Proof.* By (F1), $`z \notin L_i(x)`$ means $`x \notin L_i(z)`$, i.e. $`\kappa^{-1}(z) \notin L_i(z)`$. Apply
(D)$`[L_i(z); z]`$ for $`i = 1, 2`$. Its alternative $`\kappa^{-1}(z) \in L_i(z)`$ fails, and the second and fourth,
$`\alpha^{-1}(z) \in L_i(z)`$ and $`\sigma(z) \in L_i(z)`$, fail by (C0) at $`z`$ with $`M = L_i(z)`$. Hence
$`\rho^{-1}(x) = \rho^{-1}(\kappa^{-1} z)`$ lies in both $`L_1(z)`$ and $`L_2(z)`$, whose intersection is $`\{z\}`$ by
(F2). $`\blacksquare`$

**Lemma 3 (a composite of two OA bijections is not an OA bijection).** $`\varphi_{13} \varphi_{12}^{-1}: P_2 \to P_3`$
is not an OA bijection. Equivalently, $`\kappa \ne \mathrm{id}`$ for every OA bijection $`\varphi_{23}`$.

*Proof.* $`\kappa = \mathrm{id}`$ means $`\varphi_{23} = \varphi_{13} \varphi_{12}^{-1}`$, so the two statements are
equivalent. Suppose $`\psi = \varphi_{13} \varphi_{12}^{-1}`$ were an OA bijection. For stars $`S`$ of $`P_2`$ and $`T`$
of $`P_3`$, $`\psi(S)`$ is a PM, so $`1 = |\psi(S) \cap T| = |\varphi_{12}^{-1}(S) \cap \varphi_{13}^{-1}(T)|`$: every OA
matching of $`\varphi_{12}`$ meets every OA matching of $`\varphi_{13}`$ in exactly one cell. Now consider the six
partitions of $`P_1`$ into 4-sets

```math
\begin{array}{llll}
\text{rows}, & \text{columns}, & L_1 = \{\varphi_{12}^{-1}(\text{rows of } P_2)\}, & L_2 = \{\varphi_{12}^{-1}(\text{columns of } P_2)\}, \\
 & & L_1' = \{\varphi_{13}^{-1}(\text{rows of } P_3)\}, & L_2' = \{\varphi_{13}^{-1}(\text{columns of } P_3)\}.
\end{array}
```

Two blocks taken from two different partitions of this list always meet in exactly one cell: a row
and a column trivially; a star and an $`L`$- or $`L'`$-block because those blocks are PMs (F0); $`L_1`$ against
$`L_2`$, and $`L_1'`$ against $`L_2'`$, as in (F2); $`L`$ against $`L'`$ by the assumption. Fix a cell $`x`$. The six
blocks through $`x`$ each contain three further cells, and these six 3-sets are pairwise disjoint (a
common cell $`y \ne x`$ would put two cells in the intersection of two blocks), so $`P_1 \setminus \{x\}`$
would contain 18 cells. It contains 15. $`\blacksquare`$

Lemma 3 is the classical bound of $`n + 1 = 5`$ mutually orthogonal parallel classes on $`n^2 = 16`$
points, i.e. at most $`n - 1 = 3`$ mutually orthogonal Latin squares of order 4: the six classes would be
rows, columns and four MOLS.

## 4. Proof of the theorem

**Proposition.** If $`\alpha, \rho, \sigma`$ satisfy (A)-(D) for some permutation $`\kappa`$ of $`P_1`$, then
$`\kappa`$ is the identity.

*Proof.* Let $`x`$ be a cell, $`z = \kappa(x)`$, and suppose $`z \ne x`$. By (F2), $`z`$ lies on at most one of the
four lines of $`x`$. Three cases.

*(a) $`z`$ lies on an OA matching $`M = L_i(x)`$ of $`x`$.* Then $`z`$ is in neither $`R(x)`$ nor $`C(x)`$, so
$`\alpha(x) = z`$ by Lemma 1, i.e. $`\alpha^{-1}(z) = x \in M`$. But $`z \in M`$, so (C0) says
$`\alpha^{-1}(z) \notin M`$. Contradiction.

*(b) $`z`$ lies on a star $`S`$ of $`x`$ ($`S = R(x)`$ or $`C(x)`$).* Then $`z`$ is in neither $`L_1(x)`$ nor $`L_2(x)`$,
so $`\rho^{-1}(x) = z \in S`$ by Lemma 2. But $`x \in S`$, so (A0) says $`\rho^{-1}(x) \notin S`$. Contradiction.

*(c) $`z`$ lies on none of the four lines of $`x`$.* Lemmas 1 and 2 both apply: $`\alpha(x) = \rho^{-1}(x) = z`$.
Look at (A)$`[R(z); x]`$. Its alternative $`x \in R(z)`$ fails (since $`z \notin R(x)`$), while the
second and third, $`\rho^{-1}(x) \in R(z)`$ and $`\alpha(x) \in R(z)`$, both hold. Two alternatives out of
"exactly one". Contradiction.

So $`\kappa(x) = x`$ for every cell $`x`$. $`\blacksquare`$

*Proof of the theorem.* A solution of (A)-(D) would make $`\kappa`$ the identity by the Proposition, but
$`\kappa = \varphi_{13}^{-1} \varphi_{23} \varphi_{12}`$ is never the identity for OA bijections by Lemma 3.
$`\blacksquare`$

The theorem holds for every triple of OA bijections, not only for the orbit representatives checked in
section 6; no symmetry reduction is used.

## 5. What the proof used, and what it did not

- **All four families are needed, but only a few instances of each.** For one cell $`x`$ with
  $`z = \kappa(x) \ne x`$ the proof reads five constraints in cases (a) and (b) and nine in case (c):
  (A)$`[R(x); x]`$, (A)$`[C(x); x]`$, (B)$`[R(x); z]`$, (B)$`[C(x); z]`$ for Lemma 1; (C)$`[L_1(z); z]`$,
  (C)$`[L_2(z); z]`$, (D)$`[L_1(z); z]`$, (D)$`[L_2(z); z]`$ for Lemma 2; and (C)$`[M; z]`$, (A)$`[S; x]`$ or
  (A)$`[R(z); x]`$ for the final clash. (A0) and (C0) use the "at most one" half of the constraints,
  Lemmas 1-2 the "at least one" half. No family can be dropped altogether: there are $`\varphi_{23}`$ for
  which (A)+(B)+(C) is satisfiable (32, 26 and 32 of the 6912 $`\varphi_{23}`$ for the three $`\varphi_{13}`$
  representatives) and $`\varphi_{23}`$ for which (A)+(C)+(D) is satisfiable (8, 34 and 32 of 6912); for
  representative 0 every pair of families is satisfiable for some $`\varphi_{23}`$, and (A)+(C) for all
  of them.
- **Nothing about $`\kappa`$ except $`\kappa \ne \mathrm{id}`$.** Lemmas 1, 2 and the Proposition hold for any
  permutation $`\kappa`$ of $`P_1`$. The OA structure of $`\varphi_{23}`$ enters only via Lemma 3.
  This is sharp: with $`\kappa = \mathrm{id}`$ the four families are satisfiable (`proof120/verify_general.py`
  prints a solution), so "$`\kappa = \mathrm{id}`$" is exactly what (A)-(D) say about $`\kappa`$. In the graph,
  $`\kappa = \mathrm{id}`$ reads $`\varphi_{13} = \varphi_{23} \varphi_{12}`$: the three OA bijections would have to form
  a commuting triangle, and Lemma 3 says the composite of two OA bijections always sends some
  hub star to a set hitting a hub of the far block twice.
- **Only the incidence facts (F0)-(F2).** The proof never uses the fifth parallel class of the affine
  plane behind $`\varphi_{12}`$ (the four common transversals of $`L_1`$ and $`L_2`$, which complete rows,
  columns, $`L_1`$, $`L_2`$ to an affine plane of order 4), nor that
  $`\alpha, \rho, \sigma`$ preserve it (they do: (A0) and (C0) applied at every cell push
  $`\alpha(x), \rho(x), \sigma(x)`$ off all four lines of $`x`$, hence onto the fifth line through $`x`$), nor any
  property of the three $`\varphi_{13}`$ representatives (their planes are: the plane of $`\varphi_{12}`$ with
  $`L_1' = L_1`$, $`L_2' = L_2`$; the other plane through rows and columns; the plane of $`\varphi_{12}`$ with
  $`L_1' = L_1`$ and $`L_2'`$ the fifth class). Orthogonality of the two Latin squares behind $`\varphi_{12}`$ is
  used only as $`L_1(x) \cap L_2(x) = \{x\}`$, and inside the proof of Lemma 3.
- **Hub-freeness** ($`\alpha, \rho, \sigma`$ move every cell off its row and column) is (A0); (C0) is its
  twin for the OA matchings. Both are used only through the lines of the one cell $`x`$.

## 6. Computational checks

All scripts are in ../sat_closure/proof120/ and work from that directory with a Python that has
python-sat (CaDiCaL) installed; each prints its findings.

- `verify_proof.py {0,1,2}`: for each $`\varphi_{13}`$ representative and all 6912 $`\varphi_{23}`$, Lemma 3
  ($`\varphi_{13} \varphi_{12}^{-1}`$ is not among the 6912 OA bijections, so $`\kappa \ne \mathrm{id}`$), and for every
  cell $`x`$ moved by $`\kappa`$ the exact constraint set of its case (a), (b) or (c) is
  unsatisfiable alone (103,680 moved cells per representative: 41,472 in case (a), 41,472 in (b),
  20,736 in (c)). For the 300 lowest-indexed $`\varphi_{23}`$ also Lemma 1 (the four (A)/(B) constraints
  force $`\alpha(x) = \kappa(x)`$) and Lemma 2 (the four (C)/(D) constraints force
  $`\rho^{-1}(x) = \kappa(x)`$) are checked as stated, by adding the negated conclusion.
- `verify_general.py`: (A)-(D) with $`\kappa = \mathrm{id}`$ is satisfiable (solution printed); with $`\kappa`$ an
  arbitrary non-identity permutation it is unsatisfiable for all 120 transpositions, 500 random
  permutations and 53 random 3-cycles; $`\varphi_{13} \varphi_{12}^{-1}`$ is not an OA bijection for the fixed
  $`\varphi_{12}`$ and all 6912 $`\varphi_{13}`$ (by the $`\mathrm{Aut} \times \mathrm{Aut}`$ symmetry this is Lemma 3
  for every pair); and the satisfiable counts for pairs of families (representative 0): (A)+(B) 776,
  (A)+(D) 200, (B)+(C) 200, (B)+(D) 200, (C)+(D) 48 of 6912, (A)+(C) all 6912.
- `triples.py ABC`, `triples.py ACD`: the satisfiable counts for the two family triples quoted in
  section 5.
- Deletion-based minimal unsatisfiable subsets of (A)-(D) at one moved cell have 5 to 13 constraints,
  and the 5-constraint cores are exactly the constraint sets of case (a); the script for these cores
  and for the affine-plane classes behind $`\varphi_{12}`$ is not included.
