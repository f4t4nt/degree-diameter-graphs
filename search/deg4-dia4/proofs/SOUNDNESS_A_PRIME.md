# Soundness of the refined excess matching ((a'), (b), excess budget) used for the 108 shape lists

This is the one-page argument for the refined excess matching: what `matching_excess`
(`windows108.py` in ../sat_closure/) assumes, and why every graph of the architecture
satisfies it. It is a proof that the filter is a *necessary* condition; nothing here claims that
surviving shapes are realisable.

## Setting and assumptions (exactly those of TEMPLATE_CALCULUS.md 13.1)

$`G`$ is a simple 4-regular graph of diameter 4 containing three vertex-disjoint subdivided $`K_{4,4}`$
blocks: hubs are adjacent only to their four ports, every port is adjacent to its two hubs and to
exactly two further vertices (its *holders*, both connectors, since there is no port-port edge), and
the remaining vertices are connectors, adjacent only to ports and connectors. For a connector $`c`$,
$`x(c)`$ is its port profile and $`f_i(c) = x_i(c) + \sum_{c' \sim c} x_i(c')`$, the sum over the connector
neighbours $`c'`$ of $`c`$ (the *carried vector*). A *unit* is a complementary triple pair
($`x + x' = (2,2,2)`$); each member's only connector neighbour is its partner, so $`f = (2,2,2)`$ for both
members and a unit holds exactly two ports of every block. Every other component of the connector
graph is a *gadget*; the shape fixes the types and the connector edges of every gadget, hence every
$`f(c)`$ of every gadget connector. Nothing else is assumed: no tightness, no bound on component
potentials, no unit structure.

## The three facts

**Fact 1 (routes).** For a hub $`a`$ of block $`i`$ and a port $`q`$ of block $`j \ne i`$, every walk of length
$`\le 4`$ from $`a`$ to $`q`$ has the form $`a - p - c - q`$ or $`a - p - c - c' - q`$ with $`p`$ a port of $`a`$, $`c`$ a
holder of $`p`$ and $`c'`$ a holder of $`q`$ (section 3: from $`p`$ the walk can continue to a hub of block $`i`$,
after which only block-$`i`$ ports are reachable in the remaining two steps, or to a holder $`c`$ of $`p`$;
from $`c`$ the walk reaches $`q`$ directly, or through a connector neighbour $`c'`$ of $`c`$ that holds $`q`$; a
port after $`c`$ would have to be adjacent to $`q`$, and ports are not adjacent).

Hence, writing $`h_1, h_2`$ for the holders of $`q`$ and $`N[h] = \{h\} \cup (\text{connector neighbours of } h)`$,
the hubs of block $`i`$ within distance 4 of $`q`$ are exactly the hubs incident to the block-$`i`$ ports of
the connectors in $`N[h_1] \cup N[h_2]`$. Let $`M_i(q)`$ be the *multiset* of block-$`i`$ ports held by the
connectors of $`N[h_1]`$ and of $`N[h_2]`$, each connector counted once for every closed neighbourhood it
belongs to. Then

```math
|M_i(q)| = f_i(h_1) + f_i(h_2),
```

and diameter 4 requires the *set* underlying $`M_i(q)`$ to cover all eight hubs of block $`i`$, so it has
at least four distinct ports; in particular $`f_i(h_1) + f_i(h_2) \ge 4`$ (the inequality port lemma).

**Fact 2 (the total excess is fixed).** Put $`e(q, i) = f_i(h_1(q)) + f_i(h_2(q)) - 4 \ge 0`$. Summing over
the 16 ports $`q`$ of block $`j`$, each connector $`c`$ holding $`x_j(c)`$ ports of block $`j`$ contributes
$`x_j(c) f_i(c)`$, so $`\sum_q e(q, i) = \sum_c x_j(c) f_i(c) - 64 = \mathrm{row}_{ji} - 64`$, and summing over the
six ordered block pairs, $`\sum_{(q, i)} e(q, i) = (B(G) - 768)/2 = \Phi + D`$ by the budget identity of
13.2 ($`B = 768 + 2(\Phi - \Phi_{\min})`$, $`\Phi_{\min} = -D`$). So the total local excess of a member of
potential $`\Phi`$ is exactly $`\Phi + D`$ (0 at budget 768, 2 at 772), with no choice involved.

**Fact 3 (exact pairs).** If $`e(q, i) = 0`$ then $`|M_i(q)| = 4`$ while its underlying set has at least four
elements, so $`M_i(q)`$ has no repeated element and is a perfect matching of $`K_{4,4}(i)`$. Consequently:

- (a') if $`h_1`$ and $`h_2`$ are adjacent, every block-$`i`$ port of $`h_1`$ is in $`M_i(q)`$ twice (through
  $`N[h_1]`$, as $`h_1`$ itself, and through $`N[h_2]`$, as a neighbour of $`h_2`$), and likewise for $`h_2`$; so
  $`x_i(h_1) = x_i(h_2) = 0`$.
- (b) a common connector neighbour $`c`$ of $`h_1`$ and $`h_2`$ lies in $`N[h_1]`$ and in $`N[h_2]`$, so
  $`x_i(c) = 0`$.

Rule (a) "co-holders of an exact pair are never adjacent" (the rule of the 106 and 107 lists) is the
special case $`x_i(h_1) + x_i(h_2) > 0`$; it fails for two adjacent connectors with no block-$`i`$ port
that co-hold a block-$`j`$ port (for instance two same-block doubles $`s_j`$, or $`s_j`$ adjacent to $`d_{jk}`$,
at the block $`i`$ outside both profiles). (a') is exactly the condition that Fact 3 proves. (c) "two
connectors co-hold at most one port" is not implied by anything here and is not used (no list uses
it; see TEMPLATE_CALCULUS.md 13.3).

## What the program checks

`matching_excess(gadgets, u, excess)` with excess $`= \Phi + D`$ receives the gadget slots (one entry per
(connector, block, copy), `shapes107.analyse`), calls a slot *weak* when its connector's carried
vector differs from 2 in some foreign coordinate, and searches exhaustively for an assignment that
gives every weak slot a partner: either another gadget slot of the same block on a different
connector, admissible when $`f_i(h_1) + f_i(h_2) \ge 4`$ for both foreign $`i`$, with cost
$`\sum_i (f_i(h_1) + f_i(h_2) - 4)`$ and, at every coordinate with sum exactly 4 and both connectors in
the same gadget, (a') and (b) as stated above; or a unit slot, admissible when both foreign
coordinates are $`\ge 2`$, with cost $`f_i + f_i' - 4`$, at most $`2u`$ unit slots per block. Each slot is
used at most once and the costs must sum to at most excess.

## Why every member passes

Take a member $`G`$ whose shape is ($`u`$ units, the given gadgets) and let the true holder pairs define,
for every weak gadget slot, its partner slot (the other holder of the same port). The partner is on
a different connector (a simple graph has no double edge), in the same block, and is a gadget slot
or a unit slot. Fact 1 gives $`f_i(h_1) + f_i(h_2) \ge 4`$ at both foreign coordinates; the program's cost
of the pair is exactly $`e(q, i)`$ summed over the two foreign coordinates; Fact 3 gives (a') and (b) at
every exact coordinate (when the two holders are in different components they are neither adjacent
nor have a common neighbour, so the program correctly checks nothing there, and a unit member's only
neighbour is its partner, so the same holds for unit partners). A unit partner is possible only if
$`f_i(h) + 2 \ge 4`$ for both foreign $`i`$, which is the program's condition, and a unit holds two ports
per block, so at most $`2u`$ gadget slots per block have unit partners. Finally the true pairs are
disjoint (each slot is one port incidence), and the sum of the program's costs over the represented
pairs is a sub-sum of $`\sum_{(q, i)} e(q, i) = \Phi + D`$ (Fact 2), so it is at most excess. Hence the
true pairing restricted to the weak slots is one of the assignments the exhaustive search examines,
and the function returns True. QED.

## Consequences

1. The refined matching is sound at every budget. The 4,151 shapes of the inequality-only 772 list
   that fail it (the "inequality-only" shapes) are excluded by a necessary condition and need no SAT
   decision. The 2,325 refined shapes are the complete 772 list.
2. The 768 list itself depends on this argument: `windows108.py shapes --noexcess` applies the refined
   matching to tight shapes (excess 0) in every pipeline, and it removes 4,259 realisations at 768.
   The argument above covers excess 0 as a special case.
3. The 106, 107 and 108/776-and-above lists use (a); its soundness there rests on the separate
   argument that components of potential $`\le -10`$ do not occur at those budgets
   (TEMPLATE_CALCULUS.md 13.3).
4. Consistency checks, not proofs: all 562 S8 shapes at 772 that fail the refined matching are UNSAT by
   typed SAT anyway; the end-to-end 106 pass of the window pipeline with the refined matching keeps
   all 51 shapes of budget 768 and the 776 shape of the 106 list and removes exactly its two 772 shapes
   (both UNSAT by typed SAT), i.e. the filter removes no shape that typed SAT does not also exclude.
