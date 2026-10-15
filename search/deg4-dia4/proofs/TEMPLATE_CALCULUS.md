# Template calculus for the three-K4,4-block family

Every statement below is either proved here or checked by a named script in ../template_calculus/
or ../sat_closure/ (its outputs land under runs/ next to the script). Nothing here searches over graphs.

## 1. The family

A member is a simple graph of maximum degree 4 and diameter 4 containing three vertex-disjoint
subdivided K4,4 blocks. Block i has hubs L_i, R_i (4 + 4), ports P_i (16, one per hub pair (l, r)),
and the hub edges l - p_{lr} - r. Hubs are saturated by block edges. Every port has two free slots.
All remaining vertices are connectors; every non-block edge joins a port to a connector or two
connectors (no port-port edges; the only relaxation of this assumption is noted in section 7).

A connector c has port profile x(c) = (x_1, x_2, x_3), x_j = number of ports of block j adjacent to
c, and 4 - |x(c)| connector edges, |x| = x_1 + x_2 + x_3. Slots: sum_c x_j(c) = 32 per block. Names:
quad |x| = 4, triple 3, double 2, single 1, router 0. A double is cross ({i,j}) or same-block.

`verify.py corpus`: all 35 known perfect-104 classes are members, each with
exactly one three-block decomposition; all connectors are triples with one connector edge (a perfect
matching); class 29 is the pure (1,1,1) template, class 17 has six bent (2,1,0)-type triples, the
other 33 have four; every bent triple is matched to a complementary bent triple, so every matched
pair has exactly two ports per block.

## 2. Within-block distances and why the blocks are K4,4

Inside a block: hub-hub 2 (opposite sides) or 4 (same side); hub-port <= 3 (a hub is adjacent to one
endpoint of every edge of K4,4); port-port 2 (shared hub) or 4. So every within-block pair is within
4 using block edges alone, and no connector capacity is spent inside blocks.

Lemma. K4,4 is the only 4-regular base graph with this property. The hub-port bound needs every
vertex to be adjacent to or on every edge: for an edge uv, N[u] u N[v] = V, so b <= |N[u]| + |N[v]|
- |{u,v}| <= 8, with equality only if N[u] and N[v] meet exactly in {u,v}, i.e. the graph is
triangle-free; a triangle-free 4-regular graph on 8 vertices has 16 edges, the Mantel/Turan maximum,
hence is K4,4. Any larger or different base graph has within-block pairs at distance >= 5 that must be
repaired by connector routes, which the budget below does not pay for. This is why the record's
architecture is forced to use K4,4 and why larger blocks are a worse trade, not a free one.

## 3. Path shapes between vertex classes (complete enumerations)

Hubs are adjacent only to their ports; ports to their two hubs and their slot connectors; connectors
to ports and connectors. Enumerating a-?-?-?-b by the class of each step gives, for a in block i:

- hub a -> port q of block j != i: a-p-c-q and a-p-c-c'-q only (p a port of a, c on p's slot).
- hub a -> hub b of block j: a-p-c-q-b only.
- hub a -> connector c with x_i(c) = 0: a-p-c0-c, a-p-c0-c1-c, a-p-c0-q-c (c0 on a's slot; q a port
  shared by c0 and c). If x_i(c) > 0 the pair is automatically within 4 (a-p-b-p'-c).
- port p of block i -> port q of block j: p-c-q; p-c-c'-q; p-c-c'-c''-q; p-c-y-b-q (y, q share hub
  b); p-c-y-c'-q (c, c' share port y); p-a-p'-c-q (p, p' share hub a).

## 4. The hub-to-port identity

For hub a of block i and block j != i define S_j(a) = sum over the eight slot connectors c of a's
four ports of f_j(c), where f_j(c) = x_j(c) + sum_{c' ~ c} x_j(c'). By the first bullet of section 3,
S_j(a) is the number of (slot, path) ways a reaches ports of block j, counted with multiplicity, so
diameter 4 requires S_j(a) >= 16 (16 distinct ports). Each port of c in block i contributes c to two
hubs of block i, hence

    sum_{a in block i} S_j(a) = 2 sum_c x_i(c) f_j(c)             (per ordered block pair)
    sum_{a} sum_{j != block(a)} S_j(a) = sum_c self(c) + sum_{cc'} w(c,c')   (total)

with self(c) = 2(|x|^2 - sum_j x_j^2) and w(c,c') = 2 sum_j [(|x|-x_j) x'_j + (|x'|-x'_j) x_j].
Requirements: sum_c x_i(c) f_j(c) >= 64 for each of the six ordered pairs; total >= 768.
`verify.py` evaluates both sides on the 35 classes: both equal 768 in every class, and every hub reaches
all 16 ports of each other block, i.e. the known records are exactly tight (no hub sees any port of
another block twice along these paths).

Consequences for the pure-triple template (all x = (1,1,1), perfect matching): each pair contributes
exactly 48, so 16 pairs give 768 = the requirement, which recovers the n = 13b <= 104 cap of the
anatomy study as the statement "no slack anywhere". A triple with a non-triple partner, a quad, a
router, or a non-complementary bent pair contributes strictly less than 8 per slot, while a
cross-double adjacent to a triple or to a same-block double of the block it misses contributes more
(10 per slot); so above 104 the budget is paid for by double-rich gadgets only.

## 5. Counting relaxation over port profiles

`counting_ilp.py` treats the counts of connector profiles and of connector-edge profile
pairs as integer variables with: per-block slot sums, degree sums, the total and per-ordered-pair
budget rows, hub-to-hub coverage (sum over connectors of x_i x_j >= 16 per block pair, from the
single path shape a-p-c-q-b), and for routers only the hub-to-connector demand (24 hubs, supplied by
connector edges and 2-paths, with a 2-path consistency relaxation). This is a sound relaxation.

Results: maximum n is 127 for the family; the relaxation is loose above 105 because it cannot see
port sharing, port-to-port or connector-to-connector requirements. At n = 105 the maximum budget is
exactly 768, so every 105-vertex member is hub-to-port tight, and the only x-class shape that is
feasible is 30 triples + 3 doubles (quads, singles and routers are excluded at 105 by the per-pair
rows: `shapes105.py` lists the 119 infeasible shapes).

An earlier version of this relaxation used a hub-to-connector constraint for all connectors that
omitted the shared-port path a-p-c0-q-c; it gave n <= 120 and zero slack at every order. That
constraint was invalid and those conclusions are withdrawn; only the statements in this document
stand.

## 6. n = 105: exhaustive gadget analysis

`shapes105.py`. Call a matched pair of triples with two ports per block
a normal pair (both members have f = (2,2,2)); everything else is special. For the single feasible
shape it enumerates every labelled configuration of the three doubles (profiles, connector edges among
them) and of the triples attached to them (profile multisets per owner), 1728 skeletons, and applies:

1. the six per-pair budget rows exactly (a normal pair contributes 4 to each row);
2. per-hub compensation: a hub needs sum f_j >= 16 over its eight slots for both other blocks; a
   normal-pair triple contributes exactly 2, so the special incidences of a block, with deficits
   f - 2, must be placeable in the block's 4x4 port grid (two slots per port, distinct ports per
   connector) with all row and column sums >= 0 in both coordinates (binary feasibility program);
3. hub-to-connector with the shared-port path, unknown sharers bounded above;
4. hub-to-hub coverage.

246 configurations pass the budget rows, 162 of those fail per-hub compensation, 84 pass everything.
Up to block relabelling and ordering these are 17 gadget types:

- triangles of three doubles: {d12, d13, d23}; {s1, s2, s3} (same-block doubles in three blocks);
  {d_ij, d_ij, s_k};
- paths t - d_ij - s_k - d_ij - t with two straight stub triples;
- three same-block doubles s1, s2, s3 with bent stub triples: no connector edges (each s_i carries two
  bent stubs whose profiles put two ports in the blocks s_i misses), one edge, or a path.

So counting does not exclude 105 in the family; it reduces any 105 to "14 or 15 normal pairs plus
one of these gadgets" with hub-to-port tightness everywhere.

`surgery.py t1`: replacing one matched pair of a known 104 by the
triangle {d12, d13, d23} in all 8 orientations, for all 16 pairs of all 35 classes (4480 graphs),
never gets below 10 pairs at distance 5 (class 0 and class 6), typically 13-17; the failing pairs
are dominated by connector-double and hub-double pairs. The known records are therefore not one
local move away from 105 through this gadget, which agrees with the campaign's radius-2 closures.

## 7. The port lemma (tightness localised)

In a tight member every hub a of block i reaches every port q of block j != i exactly once, along
a-p-c-q or a-p-c-c'-q. The connectors that carry q within one connector step are
N[holders(q)] = the two slot connectors of q together with their connector neighbours; a hub of block
i reaches q once for every block-i port carried by a connector of that set. Hence:

    for every port q and every block i != block(q), the block-i ports of N[holders(q)] are exactly
    four distinct ports forming a perfect matching of K4,4(i).

A normal pair (two triples, two ports per block) carries exactly its two block-i ports when one of
its members holds q. `port_lemma.py` applies this to the 84 configurations of
section 6: for each one it enumerates how the gadget slots can share ports (with a normal connector,
or with another gadget slot of the same block) and demands count four and pairwise-distinct carried
ports for every gadget-held port and every other block. Result: every path, one-edge and no-edge
gadget dies; exactly the 21 triangle configurations survive (the three triangle types up to block
relabelling). For a triangle gadget {x, y, z} the closed neighbourhood of each member is the whole
triangle, so N[holders] contributes the triangle's two ports in each block: a triangle is a
"virtual pair" with the same reach sets as a normal pair.

## 8. Structure theorem for tight designs

Call a unit a normal pair or a virtual pair; a unit owns two ports in each block. Apply the port
lemma at a port q held by units U, U': the block-i ports of U and U' together form a perfect
matching, so each unit's two block-i ports are disjoint K4,4 edges (a 2-matching) and the two units
are complementary in block i. On the 72 two-matchings of K4,4, complementarity (disjoint, union a
perfect matching) is a disjoint union of eighteen 4-cycles K_{2,2}: the 4-cycle of the hub pairings
{l1,l2}{l3,l4} and {r1,r2}{r3,r4} consists of the two matchings on {l1,l2}x{r1,r2} and the two on
{l3,l4}x{r3,r4}, and its sibling uses the complementary eight edges. Since every unit meets other
units at four ports outside block i and all of them must be complementary to it in block i, all 16
units lie on one sibling pair of 4-cycles of block i (one L-pairing and one R-pairing per block),
eight units on each 4-cycle, two units per 2-matching. The two units with the same 2-matching share
both block-i ports: the unit graph of block i is a perfect matching of units with every edge doubled
(twins tau_i), and twins in block j must have complementary labels in block i. Equivalently the
design is three fixed-point-free involutions tau_1, tau_2, tau_3 on 16 units such that for each
block i the orbits of <tau_j, tau_k> split 8 + 8 into the two 4-cycles, their alternating bipartition
gives the sides, and tau_i preserves sides. `unit_structure.py` confirms every prediction on the 35
known classes (twins in every block, orbits of size 4).

The symmetry group of a block's labels is the group induced by the K4,4 automorphisms that preserve
the two hub pairings; it has order 64 and does not contain a single within-side swap (A1 <-> A2 alone),
because every pairing-preserving hub transposition swaps the labels of two sides at once (l1 <-> l2
swaps within A and within C). A first version of the enumeration used the full automorphism group of
the label hierarchy (order 128) for its normalisation and dedup; that merged non-isomorphic
skeletons and dropped half of the label patterns, and it was caught when the enumerated 104
realisations of the chosen representative failed to contain 26 of the known classes. Those earlier
numbers (15 skeletons, "all known classes are skeleton 13") are withdrawn.

`unit_designs2.py` redoes it with the true group: block 1 canonical, tau_2 up to the true block-1
stabiliser (order 16384, four orbit representatives), tau_3 up to the stabiliser of tau_2, every
side structure, every label pattern reduced per block by the 64-element label group, and finally
nauty certificates (pynauty) of the unit-port-hub incidence graph, whose automorphisms are exactly the
design symmetries: 1328 raw designs, **54 isomorphism classes of skeletons**
(sat_closure/data/skeletons54.json). `known_designs2.py` certifies each known class's incidence graph
directly from its own blocks and ports: the 35 known classes are splittings of three skeletons,
50 (8 classes), 51 (18) and 53 (9).

## 9. Exact completion (CaDiCaL through python-sat, `sat105.py`)

For a skeleton, each unit is realised as one of 10 splits (connector x holds the unit's first block-1
port and two of the other five, y the rest), or, in the 105 case for exactly one unit, as one of the
15 triangle gadget options (cross triangle, 8 orientations; same-block triangle; two cross doubles
plus a same-block double, 6). Option selectors with exactly-one; adjacency literals defined from the
selectors; for every vertex pair not within a block a clause over all witness walks of length <= 4 in
the potential graph (each witness the conjunction of its variable adjacency literals; walks through
fixed hub-port edges carry no literal). About 92k variables and 290k clauses per instance. Every reported solution is decoded and re-verified by independent BFS.

- 104, no gadget (`sat105.py 104 all`): exactly the three skeletons used by the known classes
  (50, 51, 53) are SAT, each with a verified 4-regular diameter-4 graph on 104 vertices; the other
  51 skeletons are UNSAT. Every diameter-4 member of the family at 104 is therefore a
  split of one of those three skeletons.
- 105, triangle gadgets (`sat105.py 105 all`): all 54 skeletons UNSAT (the four used by the
  known classes and their neighbours, 50 to 53, are the slowest; the rest are immediate).

Validation (`validate.py`): for each of the 35 known classes the encoder built
on that graph's own unit design, with every option fixed to the split the graph actually uses, is
SAT (35/35), so no clause excludes a genuine diameter-4 member; a gadget control in 105 mode with
only hub/port pairs encoded places a cross triangle and yields a well-formed 4-regular 105-vertex
graph with all hub/port pairs within 4, exercising the gadget options on a satisfiable instance.
The enumeration of all 104 realisations of the three live skeletons with nauty classification
(enumerate104.py, classify104.py) recovers every known class and is the independent end-to-end check
of the skeleton theory: it returns 36 classes, 8 + 18 + 10, the 35 known ones and one new one
(class 35 of ../classes/).

## 10. Theorem and its proof chain

**No simple graph with maximum degree 4, diameter 4 and 105 vertices contains three vertex-disjoint
subdivided K4,4 blocks with no edge between ports.**

1. Sections 3 and 4: the hub-to-port identity and the per-ordered-pair rows are necessary.
2. Section 5: the counting relaxation at 105 has maximum budget exactly 768, so a 105 member is
   tight, and the only x-class shape is 30 triples + 3 doubles.
3. Section 6: exhaustive enumeration of special-connector configurations under exact budget rows,
   per-hub compensation, hub-to-connector and hub-to-hub checks leaves 84 configurations.
4. Section 7: the port lemma leaves the 21 triangle configurations; triangles are virtual pairs.
5. Section 8: a tight design of 16 units is one of 54 skeletons (complete orderly enumeration with
   the true symmetry group and nauty certificates).
6. Section 9: for each skeleton, the completion over all splits and all triangle options is UNSAT.

Every step is either a proof in this document or a terminating exact computation whose inputs and
outputs are reproducible from the scripts named. The encoder is validated on all 35 known graphs with fixed splits, on the
three live skeletons at 104, and by the recovery of the known classes from the full enumeration.

## 11. The 120-vertex all-cross-double shape

The counting relaxation's extreme structure is n = 120 with 48 cross-doubles, 16 of each type
{1,2}, {1,3}, {2,3}, and no other connectors. Its budget is exactly 768 (self 4 each, every
connector edge between different types worth 12, same-type edges worth 8 and therefore excluded), so
any such member is tight. Consequences, all forced:

- hub-to-hub exactness: the 16 {i,j}-doubles cover the 64 hub pairs of blocks i, j exactly once, so
  every port holds exactly one double of each of its two types, and the {i,j}-doubles define an
  "OA bijection" phi_ij: P_i -> P_j mapping every hub star to a perfect matching (a pair of
  orthogonal Latin squares of order 4). There are 6912 labelled OA bijections between two labelled
  K4,4 blocks, a single orbit under Aut(K4,4) x Aut(K4,4) (`design120.oa_bijections`); with phi_12
  fixed, phi_13 has three orbit representatives.
- the port lemma at a port of block j held by d_ij and d_jk gives 5 - e carried block-i ports, where e
  is the number of {j,k}-neighbours of d_ij, so e = 1: every double has exactly one neighbour of each
  other type, and the connector graph is a union of cycles with type period 3.
- if all cycles are triangles, each triangle is a virtual pair, so the design is one of the 54
  skeletons with every unit a triangle: `sat105.py 120 all` is UNSAT for
  all 54 immediately.
- in general (`design120.py`): phi_12, phi_13 fixed up to symmetry, phi_23 and
  the three perfect matchings between type classes as variables, and the twelve families of
  hub-to-port exactly-one constraints (72k variables, 207k clauses): UNSAT for all three phi_13
  representatives, instantly.
- independent cross-check (`design120_reduced.py`): transporting
  everything to the ports of block 1 gives three permutations alpha, rho, sigma and a fixed
  kappa = phi_13^-1 phi_23 phi_12; hub exactness for the hubs of blocks 1 and 2 alone reads, for
  every star S and every OA perfect matching M,
      S | rho S | alpha^-1 S | sigma S = P_1,      S | alpha S | kappa rho S | kappa sigma S = P_1,
      M | alpha M | rho^-1 M | sigma^-1 M = P_1,   kappa M | alpha M | kappa rho M | sigma^-1 M = P_1.
  Over all 3 x 6912 choices of (phi_13, phi_23) these four families are UNSAT (20,736 instances),
  while each single family and the pairs (A,B), (C,D) are satisfiable, so the contradiction needs
  the hubs of two blocks together.

**Theorem.** No 120-vertex graph of maximum degree 4 and diameter 4 consists of three subdivided
K4,4 blocks and 48 cross-doubles. Together with section 10 the two budget-tight shapes of the family
(the record's own shape at 105 and the all-double extreme at 120) are closed.

**Human proof (PROOF120.md).** The four families force kappa = id: for a cell x moved by kappa,
z = kappa(x) lies on at most one of the four lines of x (row, column, two OA matchings); (A)+(B)
force alpha(x) = z whenever z is off the row and column of x, (C)+(D) force rho^-1(x) = z whenever z
is off the OA matchings, and each of the three cases then contradicts one exactly-one constraint
(at most nine constraints used in total). But kappa = id means phi_13 = phi_23 phi_12, and the
composite of two OA bijections is never an OA bijection: it would give six mutually orthogonal
partitions of the 16 ports into 4-sets, i.e. four MOLS of order 4, where at most three exist. So the
120 design fails for exactly one reason: the three orthogonal-array bijections between the blocks
would have to form a commuting triangle, and orthogonal Latin squares of order 4 cannot do that.
Every lemma was re-verified computationally for all 3 x 6912 instances (sat_closure/proof120/).

## 11b. Four octahedron blocks at 108 (sat_closure/oct4/)

Free blocks (every internal pair within 4) are exactly the subdivided K5, octahedron and K4,4
(blocks_relax.py). Four octahedra give 72 block vertices and 96 slots, the same as three K4,4, and
the hub-reach relaxation is infeasible below 108 and has zero slack at 108, so an Oct x4 member on
108 vertices would again be tight. A potential-function argument (T(c) = |x||F| - x.F summed over
connectors, redistributed over connector edges; finite check over all 126 connector types and 1596
type pairs, ILP cross-check) shows the only tight structures are 12 "unit" pairs of triples with
profile sum a permutation of (2,2,1,1), alternating cycles of complementary cross doubles, and
same-block doubles between two straight triples. The holder lemma (the two holders of a port carry
exactly three pairwise disjoint octahedron edges per foreign block, a perfect matching of K_{2,2,2})
gives carried totals {4, 5} at every port, which kills same-block doubles and makes every port
have one heavy and one light holder. Finally, a unit heavy in blocks j and k has two block-j
ports; at each, the light co-holder must carry the unique edge completing the unit's block-k
2-matching to a perfect matching (an octahedron 2-matching has at most one completion), so both
light co-holders are the single light holder of that edge, hence the same entity, which then holds
two block-j ports although light entities hold one. **No Oct x4 graph on 108 vertices exists.**
Typed SAT over the 11 possible connector skeletons confirms: all UNSAT (three encodings), encoder
validated on the known 104.

## 12. What remains open

- 106 to 119, and the slack shapes at 120 (routers, singles, quads appear in the relaxation's optima
  there with budget up to 800). The relaxation has slack from 106 on (780 at 106, growing), so
  tightness is no longer forced and the gadget zoo is larger; the shape enumeration needs
  sharing-aware hub-to-connector and port-to-port rows (section 3), and the port lemma must be used
  in its inequality form (at least four carried ports covering all eight hubs).
- Port-port edges. The family excludes them; they remove slots, so they lower n for a fixed connector
  mix, but a complete theorem should treat them. The path enumerations of section 3 gain extra shapes
  through such edges and the identity of section 4 needs the corresponding terms.

## 13. Closure of the three-K4,4 architecture at 105-107 and the slack sub-families at 108

(This section collects the proof chains of the slack closures and updates the first item of
section 12: 106 and 107 are closed, and so is 108 for every budget above 772. Scripts are in
../sat_closure/ unless noted; the shape lists are in ../sat_closure/data/.)

**13.1 Architecture.** A member is a simple 4-regular graph G of diameter 4 on n = 72 + N vertices
that contains three vertex-disjoint subdivided K4,4 blocks: block i has hubs L_i, R_i (4 + 4),
ports p_lr (16) and the edges l - p_lr - r; hubs are saturated by block edges, every port has two
further neighbours (its slots, 96 in all), there is no edge between two ports, and the remaining N
vertices are the connectors. A connector c has profile x(c) in Z^3 (x_j = block-j ports adjacent
to c), |x| <= 4, and 4 - |x| connector neighbours: triples (|x| = 3; pure (1,1,1) or bent), cross
doubles d_ij = e_i + e_j, same-block doubles s_k = 2 e_k, singles, routers (|x| = 0) and isolated
quads. Nothing else is assumed (no tightness, no unit structure). 4-regularity with every slot
filled is the family definition of section 1 and is used by every count below, including the
105 chain of section 10 (counting_ilp.py fixes the slot sums to 32 and the connector degrees to
4 - |x|), so the theorems of this section are stated for 4-regular graphs.

**13.2 Budget identity and potential.** For a hub a of block i and a block j != i let S_j(a) be the
sum of f_j(c) = x_j(c) + sum_{c' ~ c} x_j(c') over the eight slot connectors c of a's four ports.
By section 3 the routes from a to the ports of block j are a-p-c-q and a-p-c-c'-q, so diameter 4
forces S_j(a) >= 16, and summing over the 48 (hub, foreign block) pairs gives the budget
B(G) = sum_{a,j} S_j(a) = 2 sum_c T(c) >= 768, T(c) = |x||F| - x.F, F(c) = x(c) + sum_{c'~c} x(c'),
with the six rows sum_c x_i(c) f_j(c) >= 64 (section 4). Put phi(c) = T(c) - 3|x(c)| - 3 and
Phi = sum_c phi(c). Since sum_c |x(c)| = 96, B = 768 + 2 (Phi - Phi_min) with Phi_min = 96 - 3N,
i.e. Phi_min = -3, -6, -9, -12 at n = 105, 106, 107, 108, and diameter 4 forces Phi >= Phi_min.
Phi distributes over connector edges: E(c,c') = 2 v(c,c') + node(c) + node(c'), v = |x||x'| - x.x',
node(c) = (self(c) - 3|x| - 3) / (4 - |x|), self = |x|^2 - x.x, and Phi = sum_edges E + sum_quads
(self - 15). The finite table `potential3.py table` (template_calculus) gives
E <= 0 for all 210 pairs of the 20 non-isolated types, E = 0 exactly for complementary triple
pairs (x + x' = (2,2,2)) and for d_ij adjacent to s_k with {i,j,k} = {1,2,3}, and self - 15 in
{-5, -7, -9, -15} for quads (re-derived with exact rationals in `check_potential_table.py`). Hence every connected typed component of the connector graph has
potential <= 0, a member of potential Phi has budget 768 + 2 (Phi - Phi_min), no component of
potential < Phi_min occurs, and the slot deficit sum_c (3 - |x(c)|) equals D = 3N - 96 (3, 6, 9, 12),
quads counting -1. Units are the complementary triple pairs (potential 0); the other potential-0
components are the alternating cycles (d_ij s_k)^m. Everything else is a gadget.

**13.3 Enumeration: components, then shapes.** Components (`components107.py`): every
connected typed connector graph with potential >= Phi_min and deficit <= D + floor(D/5) (a component
of deficit d > D forces d - D quads, each costing >= 5; this covers the quad case), grown as cores of
doubles/singles/routers with triple leaves and deduplicated by the coloured certificate (every
colour with its multiplicity in front of the pynauty certificate; the flat certificate of the
first enumerator, `potential3.py components`, merged 432 components at 106 and is superseded;
networkx isomorphism self-check over all 19,604 components at 106). Counts
(unit types included): 12,566 at 106 with the deficit-indexed bound (19,604 with the flat bound
-6, deficit <= 7), 947,141 at 107, 456,130 (Phi >= -6) and 3,387,019 (Phi >= -8) at 108. A shape
is u units plus a multiset of gadgets with total deficit D, total
potential >= Phi_min, gadget slots 32 - 2u per block, the six rows (a unit contributes 4 to each),
the port lemma in inequality form (section 7: at a port q of block j with
holders h1, h2 and every i != j, f_i(h1) + f_i(h2) >= 4; unit members have f = (2,2,2), so this is a
matching problem on the gadget slots), and, for tight shapes only (B = 768), its equality form:
exact complements at every slot, (b) a common connector neighbour of two co-holders holds no
foreign port, (a) co-holders are not adjacent (sound while every component has potential > -10,
true in all lists here). The former condition (c) "two connectors co-hold at most one port" is not
implied by exactness and was removed (`recheck_tight.py`). Shapes are
taken up to the S_3 of block relabellings, under which the SAT models are symmetric. At 108 the
class-level pipeline `shapes108.py` (gadgets grouped by their multiset of (profile, carried
vector), disk-backed store) reproduces the sub-family lists 106/Phi >= -4, 107/Phi >= -7 and
108/Phi >= -6 of the older pipeline exactly. No surviving shape contains a quad or a router; one
shape at 107 and one at 108/776 contain a single.

| n | Phi_min | components | shapes by budget (classes up to block relabelling) | total |
|---|---|---|---|---|
| 105 | -3 | (sections 5-7) | 768 only: 1 x-class shape (30 triples + 3 doubles), 84 configurations, 21 after the port lemma | 3 triangle types |
| 106 | -6 | 12,566 | 776: 1, 772: 2, 768: 51 (183 labelled shapes; 11-14 units) | 54 |
| 107 | -9 | 947,141 | 780: 6, 776: 51, 772: 99, 768: 573 (9-13 units) | 729 |
| 108 | -12 | 456,130 (>= -6), 3,387,019 (>= -8) | 792: 1, 784: 15, 780: 374, 776: 2,938 (8-12 units); 772 and 768 not enumerated | 3,328 |

(data/shapes106.json, data/shapes107.json and data/shapes108_776plus.json.)

**13.4 Deciding a shape.** Typed SAT (`typed_slack_sat.py`, CaDiCaL 1.5.3 through
python-sat): variables A[c][q] "connector c holds port q" with the gadget profiles and connector
edges fixed and the unit types free (two adjacent triples with two ports per block in total), two
holders per port; derived literals U (c holds a port of hub a), Z (within 2), S (shared port),
R3, R4; one clause per vertex pair that is not automatically within 4, over all walk shapes of
section 3 (hub-port, hub-hub, hub-connector, port-port, port-connector, connector-connector);
redundant clauses proved in sections 4 and 7 (inequality port lemma, exactness at holder pairs
with f-sum 4, unit members never share a port) and tight mode for B = 768; symmetry breaking
(first gadget connector of each block at port (0,0), second-port orbit, ordered units). 250-365k
variables, 1.0-1.45M clauses; caps are enforced with conflict budgets; every SAT answer is decoded
to an edge list and BFS-checked independently (none occurred). Skeleton method (`sat105.py`,
`sat106.py`): a tight shape all of whose components are units or triangles of doubles
(virtual pairs, section 7) is one of the 54 skeletons of section 8, and the completion over the
10 splits per unit and 15 triangle options per gadget unit is a second, independent encoding.
Validation on the 36 known 104s (the 35 campaign classes and class 35): the
sat105.py encoder with the splits fixed to each graph's own is SAT on all of them (validate.py),
the 104 skeleton run is SAT exactly on the three skeletons 50, 51, 53 used
by the known graphs and the full enumeration of their realisations returns all 36 classes
(enumerate104.py, classify104.py); `typed_slack_sat.py validate` on class 35: the graph's own port
assignment is SAT under assumptions, SAT breaks after one port swap, the free model is SAT with a
decoded BFS-verified 4-regular diameter-4 graph and stays SAT with one unit fixed as a gadget;
negative control (`typed_slack_sat.py control`) 15 units + cross triangle at 105 UNSAT, agreeing
with sat105.py; six 106 shapes UNSAT with every redundant clause and the tight mode switched off.

**13.5 Theorem.** Let G be a 4-regular graph of diameter 4 on n vertices that contains three
vertex-disjoint subdivided K4,4 blocks with saturated hubs and no port-port edge. Then n is not
105, 106 or 107, and if n = 108 then B(G) <= 772, i.e. Phi(G) is -10 or -12. Proof chain:

1. n = 105 (sections 5-10): the counting relaxation forces B = 768; 84 special-connector
   configurations, 21 survive the port lemma (three triangle types); the 54 skeletons are UNSAT
   (`sat105.py 105 all`); the typed SAT agrees on the cross-triangle shape (UNSAT).
2. n = 106: 54 shape classes (data/shapes106.json). The 92 labelled shapes of the list obtained
   with (c) (27 classes): UNSAT; the 27 further classes found without (c): UNSAT; the 15
   two-triangle labelled shapes are UNSAT a second time by `sat106.py 106 all`.
3. n = 107: 729 classes (data/shapes107.json). 728 UNSAT by typed SAT (302 non-unit-like and 26
   unit-like from the first list, 400 from the list without (c)), and the three-cross-triangle
   shape (tight, all components virtual pairs; the typed SAT did not decide it within its cap)
   UNSAT by `sat106.py 107 all` (54 skeletons).
4. n = 108, B >= 780 (Phi >= -6): 390 shapes UNSAT (the first 390 of data/shapes108_776plus.json).
5. n = 108, B = 776 (Phi = -8): 2,938 shapes UNSAT (the rest of data/shapes108_776plus.json); the
   twelve that needed a longer cap all have 12 units: ten 12-cycles of doubles, a cross triangle +
   9-cycle, a cross triangle + {d13, d13, s2} triangle + 6-cycle.

Every step is a proof in this document or a terminating exact computation; the encoders are the
validated ones of 13.4, imported unchanged by the driver (`typed_sat_shapes.py`).

**13.6 Open.** (i) n = 108 with B = 772 or 768 (Phi = -10, -12): the component families have
roughly 20M and 100M+ members, out of reach of the component pipeline; the equality form of the
port lemma applies at all but at most two (port, block) pairs (772) or everywhere (768), and the
intended route is an enumeration of connector windows (a connector with its neighbours' types)
under the Eulerian chaining condition and the exact slot-complement balance, realised as connector
graphs only afterwards (13.7). (ii) n >= 109 (Phi_min <= -15), where component
counts grow about 50x per three units of potential. (iii) Mixed architectures Oct + 2 K4,4 and
K5 + 2 K4,4: the counting slack is 48-58, there is no potential cut and hence no finite shape
list; a validated typed encoder for mixed block lists exists (not included here), and on the
restricted family R(105) (631 shapes: 3 units, one cycle per class, two leaf paths) it decided
31 shapes UNSAT with 25 undecided at short caps and no SAT answer, 575 not run; the hub-only
staged model is slower than the full model there, and closing R(105) alone would need far more
solver time than any of the closed cases. (iv) Port-port edges,
which the architecture excludes: they add route shapes to section 3 and terms to the identity of
section 4, and a complete theorem for "contains three subdivided K4,4 blocks" must treat them.

**13.7 Closure at n = 108, budgets 772 and 768.** The component
pipeline of 13.3 does not reach Phi = -10, -12 (tens of millions of components), so the
closure re-expressed a shape as a multiset of *windows* (a connector type with the multiset of its neighbours'
types): -Phi is a bounded knapsack over window costs, the Eulerian chaining condition and the slot
complements are linear balance conditions, and a window multiset is realised as typed graphs only
afterwards (windows108.py, validated by reproducing the 54 shapes of 106 at every budget with two
independent combine stages). A C program then enumerated all 7,702,003 chaining-closed window multisets
("clusters") of cost <= 12 with at most four unprovidable slots (clusters108.c; equal to the
Python growth where both apply, 3,150 = 3,150 at 106 and 143,889 = 143,889 at 108 with size <= 8, and
complete where the Python growth's seed filter was not), proved that every non-root cluster of a
shape is cheap (cost <= 6; 361,164 clusters) and combined clusters over any number of components by
a class-level group combine (combine_classes.c; pair combine C = Python, group combine contains the
pair combine, class combine = group combine on the size-8 clusters, and the whole chain reproduces
the 106 lists end to end). Realisation, canonical form under block relabelling, the inequality port
lemma and the refined excess matching give the complete lists: **11,163 shapes at 768** (6-13 units)
and **2,325 at 772** (8-12 units; 7,056 under the inequality matching alone)
(data/shapes108_768.json.gz and data/shapes108_772.json.gz).
The refined excess matching replaces condition (a) of 13.3 and is proved sound for every member at
every budget in SOUNDNESS_A_PRIME.md: the hubs of block i within 4 of a port q are those on
the block-i ports of the closed neighbourhoods of q's two holders, a multiset of size f_i(h1) + f_i(h2)
that must cover all eight hubs; the local excesses f_i(h1) + f_i(h2) - 4 sum to exactly Phi + D =
(B - 768) / 2; and at an exact pair the multiset has no repetition, so (a') adjacent co-holders hold
no block-i port and (b) a common neighbour holds none ((a) is the special case x_i(h1) + x_i(h2) > 0
and fails for adjacent co-holders without a block-i port, which first occur at potential -10).
Typed SAT (13.4 encoder unchanged, tight mode at 768; per-instance caps enforced in conflict-budget
chunks): **all 11,163 shapes at 768 and all 2,325 shapes at 772 are UNSAT**; no instance timed out
at its final cap and no SAT answer occurred (two all-triangle shapes also by the skeleton method,
`sat106.py 108 all`; the ten hardest shapes are tight 12-unit shapes whose gadgets are cycles of
doubles with type period 3, and the hardest, two cross triangles plus a cross 6-cycle, was decided
a second time by a symmetry-reduced case split, `cube_sat108e.py`; the verdict files of the batches
were merged per shape by `merge108e.py`). Finally B is always a multiple of 4: T(c) = self(c) + sum_{c' ~ c} v(c, c') with
self(c) = 2 sum_{i<j} x_i x_j even, so sum_c T(c) is even and B = 2 sum_c T(c) == 0 mod 4; at 108
the budgets are therefore 768, 772 and the >= 776 of 13.5.

**13.5, updated theorem.** Let G be a 4-regular graph of diameter 4 on n vertices that contains
three vertex-disjoint subdivided K4,4 blocks with saturated hubs and no port-port edge. Then n is
not 105, 106, 107 or **108**. (Items 1-5 as before; item 6: n = 108, B = 772 and 768, by 13.7.)

**13.6, updated.** Item (i) is closed. Open: (ii) n >= 109 (Phi_min <= -15; the window-cluster
enumeration grows by a factor of roughly 50 per three units of potential, a few 1e8 clusters at 109),
(iii) mixed architectures, (iv) port-port edges, as before.
