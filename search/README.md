# Search Methods

This directory holds the search code behind the graphs in `graphs/`, one subdirectory per method.

| directory | method | cells it produces |
|---|---|---|
| `cayley_metacyclic_anneal/` | simulated annealing of Cayley connection sets on $`\mathbb{Z}_N \rtimes_r \mathbb{Z}_M`$ with level patterns, targeted moves, one movable involution for odd degree, exact completion sweeps; a two-plane variant for $`n`$ up to $`3 \times 10^{10}`$ | [(7,9)](../graphs/deg7/dia9), [(7,10)](../graphs/deg7/dia10), [(8,7)](../graphs/deg8/dia7), [(8,8)](../graphs/deg8/dia8), [(8,9)](../graphs/deg8/dia9), [(8,10)](../graphs/deg8/dia10), [(9,7)](../graphs/deg9/dia7), [(9,8)](../graphs/deg9/dia8), [(9,9)](../graphs/deg9/dia9), [(9,10)](../graphs/deg9/dia10), [(10,7)](../graphs/deg10/dia7), [(10,8)](../graphs/deg10/dia8), [(10,9)](../graphs/deg10/dia9), [(10,10)](../graphs/deg10/dia10), [(11,7)](../graphs/deg11/dia7), [(11,8)](../graphs/deg11/dia8), [(11,9)](../graphs/deg11/dia9), [(11,10)](../graphs/deg11/dia10), [(12,7)](../graphs/deg12/dia7), [(12,8)](../graphs/deg12/dia8), [(12,9)](../graphs/deg12/dia9), [(12,10)](../graphs/deg12/dia10), [(13,7)](../graphs/deg13/dia7), [(13,8)](../graphs/deg13/dia8), [(13,9)](../graphs/deg13/dia9), [(13,10)](../graphs/deg13/dia10), [(14,8)](../graphs/deg14/dia8), [(14,9)](../graphs/deg14/dia9), [(14,10)](../graphs/deg14/dia10), [(15,7)](../graphs/deg15/dia7), [(15,8)](../graphs/deg15/dia8), [(15,9)](../graphs/deg15/dia9), [(15,10)](../graphs/deg15/dia10), [(16,7)](../graphs/deg16/dia7), [(16,8)](../graphs/deg16/dia8), [(16,9)](../graphs/deg16/dia9), [(16,10)](../graphs/deg16/dia10) |
| `cayley_exhaustive/` | complete sweeps of Cayley families up to group automorphisms: the one-involution odd-degree family, level-pattern sweeps of 3-pair sets, fully reduced single-group lists, mode F with the automorphism quotient | [(6,8)](../graphs/deg6/dia8), [(6,10)](../graphs/deg6/dia10), [(7,6)](../graphs/deg7/dia6), [(7,8)](../graphs/deg7/dia8), [(9,5)](../graphs/deg9/dia5) |
| `lift_two_vertex/` | reduced sweeps of 2-vertex voltage lifts (degree-4 shapes, degree-5 shape 103 with $`S_3`$ reduction and sub-blocks, the R1 pattern filter) over metacyclic groups with a bit-parallel evaluator | [(4,10)](../graphs/deg4/dia10), [(5,8)](../graphs/deg5/dia8), [(5,9)](../graphs/deg5/dia9), [(5,10)](../graphs/deg5/dia10) |
| `lift_three_vertex/` | the same for 3-vertex base shapes | [(4,9)](../graphs/deg4/dia9) |
| `lift_nonmetacyclic_groups/` | shape-103 lifts and one-involution Cayley graphs over groups $`K \rtimes C_M`$ built from blocks (not necessarily metacyclic), chosen by the conjugacy-class-count rule | [(5,7)](../graphs/deg5/dia7) |
| `route_chart_lifts/` | route-chart lifts: a controller (a Schreier-base lift over a small group or a Cayley graph) with a fibre $`\mathbb{F}_3^4`$ over every vertex and per-symbol charts $`(M_t, d_t)`$, swept with the controllability word filter, chart/controller co-design by annealing, mixed alphabets with central plain symbols, and 27-point cells with a per-edge type SAT at diameter 3 | [(12,4)](../graphs/deg12/dia4), [(13,4)](../graphs/deg13/dia4), [(14,4)](../graphs/deg14/dia4), [(15,3)](../graphs/deg15/dia3), [(15,4)](../graphs/deg15/dia4) |
| `group_enumeration/` | enumeration and screening of the candidate groups: metacyclic presentations, the $`k(G)`$ closed form, random-set softness screens | supporting (feeds the sweeps above) |
| `hexagon_insertion/` | subdivided flags of the generalized hexagon $`H(q)`$ plus an annealed graph on the new vertices; flag sets as group orbits for $`q = 3`$ | [(4,6)](../graphs/deg4/dia6), [(5,6)](../graphs/deg5/dia6) |
| `insert_grow/` | absorbed-deletion insertion into an existing graph: a new vertex takes edges on which no distance-$`D`$ pair depends, with greedy coverage of the vertex set by its neighbours' balls, exact CaDiCaL completion of the dependent pairs, and a selection model in which the solver chooses the deletion set and the neighbourhood together | [(11,4)](../graphs/deg11/dia4) |
| `cubic_lifts/` | exhaustive voltage lifts of small cubic base graphs over all groups of a given order | [(3,7)](../graphs/deg3/dia7) |
| `deg4-dia4/` | the (4,4) machinery: the lift-and-repair search for 104 vertices and the structural classification of its 36 classes | [(4,4)](../graphs/deg4/dia4) |

Conventions shared by all Cayley and lift code: the metacyclic group is
$`\mathbb{Z}_N \rtimes_r \mathbb{Z}_M = \langle x, y \mid x^N = y^M = 1,\ y x y^{-1} = x^r \rangle`$; the element
$`x^i y^j`$ has index $`i + N j`$; the product is
$`(i_1, j_1)(i_2, j_2) = (i_1 + r^{j_1} i_2 \bmod N,\ j_1 + j_2 \bmod M)`$; the Cayley graph has edges $`g \sim g s`$
for $`s`$ in the connection set and its inverses; a voltage lift joins $`(u, g)`$ to $`(v, g a)`$ for a dart
$`u \to v`$ with voltage $`a`$, and the lift vertex $`(v, g)`$ has index $`v\,|G| + \mathrm{index}(g)`$.

## Glossary

- defect: the number of vertices farther than $`D`$ from the root (a Cayley graph is vertex-transitive, so one BFS
  decides the diameter; a lift needs one per base vertex). Defect 0 means diameter at most $`D`$.
- level, level pattern: the level of $`x^i y^j`$ is its $`y`$-exponent $`j`$; the level pattern of a connection set or
  voltage assignment is the list of the generators' levels. The search fixes a pattern and then moves only the
  $`x`$-coordinates, because the pattern alone determines how many words of length at most $`D`$ land on each
  layer.
- level-pattern block: all assignments with a given level pattern (the $`x`$-coordinates free); a complete sweep of a
  block evaluates every $`x`$-assignment of that pattern once up to the automorphisms that preserve the pattern.
- coincidence-free pattern: a level pattern in which no two reduced words of length at most 2 land on the same vertex
  of the quotient cycle $`\mathbb{Z}_M`$ (beyond the forced commutations), so no vertices are wasted at radius 2; the
  pattern optimiser searches for such patterns with the smallest yardstick.
- yardstick: the Poisson estimate of the defect of a random assignment with a given level pattern, the sum over layers
  of $`L \exp(-w/L)`$ where $`w`$ of the words of length at most $`D`$ land in a layer of size $`L`$; the uniform
  yardstick is its value for a pattern with no structure. The ratio of the two says how much a pattern helps.
- shape 103: the 2-vertex base graph of degree 5 with one loop at each vertex and three parallel links between the two
  vertices, one link carrying the trivial voltage; the lift has $`2|G|`$ vertices. Shapes of degree 4 are named B2a (a
  loop at each vertex, two links) and haar (four links).
- sub-block and $`S_3`$ reduction: a sub-block $`(i, j)`$ of the shape-103 family fixes the levels of the two loops to
  $`i`$ and $`j`$; the $`S_3`$ reduction uses the three re-gaugings of base vertex 1 that permute the three links (an
  isomorphism of the lift) to keep one of the six equivalent link pairs, on top of the group automorphisms and loop
  inversions.
- unit-$`x`$: a connection set whose generators all have $`x`$-coordinates coprime to $`N`$ (or, in the weaker form,
  coprime to the number of involutions $`I`$), which keeps every generator from commuting with the involution.
- mode F: the annealer's complete two-slot sweep, every admissible $`x`$-value of one generator times every admissible
  $`x`$-value of another, with the level pattern, the involution and the remaining generators fixed; started from any
  state with the start values included it is a complete sweep of that family, and the automorphism quotient (`--quot`)
  evaluates one state per orbit of the automorphisms fixing the involution and the fixed slots.
- frozen state: an annealing state that the annealer cannot improve (typically defect 2 with the two far vertices
  forming a pair); single-slot or $`x`$-only completion is the exact evaluation of every change of one generator (all
  values, or only its $`x`$-coordinate) of such a state, so a diameter-$`D`$ neighbour cannot be missed.
- movable involution: in the odd-degree family the involution is $`t = x^{t_x} y^{M/2}`$; the $`x`$-coordinate $`t_x`$
  is itself a search variable (a $`t`$-move), aimed at an unreached involution of the current state when there is one.
- softness screen, lower quartile: random connection sets of a family are evaluated and the family is ranked by the mean
  defect of the best quarter of the sample (the lower-quartile mean, often relative to the uniform yardstick); softer
  families are more likely to contain a defect-0 state.
- small-complement rule: among groups of a given order, prefer those whose complement order $`M`$ is small relative to
  the cube root of the order, about $`M \approx 0.4\, n^{1/3}`$; such groups have few involutions (none for odd order)
  and a large commutator subgroup.
- class-count rule: prefer the group with the fewest conjugacy classes $`k(G)`$ in its window (the lowest commuting
  probability $`k(G)/|G|`$), which selects the (5,7) group.
- two-plane bitmap annealing: the same annealing with a BFS storing only two bits per vertex (visited, and the parity
  of the level at which it was reached), so graphs of up to about $`3 \times 10^{10}`$ vertices fit in 8 GB.
- hexagon insertion: subdividing $`m`$ flags (incidence edges) of the incidence graph of the generalized hexagon
  $`H(q)`$ with new vertices joined by a graph $`F`$ so the result is $`(q+1)`$-regular and keeps diameter 6.
