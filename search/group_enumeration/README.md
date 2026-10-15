# Group Enumeration

I enumerate and screen the metacyclic groups a window of orders offers before any sweep is started: every presentation
$`G(N, M, r, t) = \langle x, y \mid x^N = 1,\ y^M = x^t,\ y x y^{-1} = x^r \rangle`$ of a given order range up to the
known equivalences, with the size of the reduction group, an isomorphism fingerprint and a rigorous lift bound, plus the
conjugacy-class count $`k(G)`$ in closed form (the commuting probability $`k(G)/|G|`$ is the selection rule for the
(5,7) group) and a random-set softness screen of whole shape-103 families.

Files:

- `groups_meta.py`: all split and non-split metacyclic groups with $`|G| \in [\mathrm{lo}, \mathrm{hi}]`$ (odd, even or
  all orders), one presentation per proven equivalence, columns `n N M r t c Mc ord |A| fam_B2a fam_haar ninv
  fingerprint shapes`, with the lift bound $`2 M c \le 2 L_3(D)`$ that rules a group out for the degree-4 shapes.
- `classcount.py`: $`k(G) = M \sum_{e \mid N} \varphi(e) / \mathrm{ord}_e(r)^2`$ for $`\mathbb{Z}_N \rtimes_r \mathbb{Z}_M`$,
  the per-layer class counts, and `test` (the formula against brute-force orbit counting on small groups).
- `screen_groups.py`: for each group of a list, the defect histogram of random shape-103 candidates with loops in unit
  layers (compiles `../lift_two_vertex/bitlift.c` on demand; the lower quartile of the defect distribution is the
  softness statistic used to rank groups).

To run it yourself (seconds each):

```
pip install -r requirements.txt
python3 classcount.py 703 18 85
# expected first line: G = Z_703 x|_85 Z_18, |G| = 12654, ord r = 18, k(G) = 76, cp = 0.00600601
python3 classcount.py test
# expected last line: 46 of 46 brute-force cases agree (formula = Burnside = brute-force orbits = commuting pairs / |G|); closed forms (prime N, dihedral, abelian) OK
python3 groups_meta.py 9600 9700 10 odd | head -3
# expected first line: # odd-order metacyclic groups G(N, M, r, t) with |G| in [9600, 9700], D = 10, lift bound M c <= L3(D) = 1561: 85 presentations, 77 pass the bound, 46 distinct fingerprints; ...
```

The group of the (7,6) graph above, $`\mathbb{Z}_{703} \rtimes_{85} \mathbb{Z}_{18}`$, has the lowest $`k(G)`$ of its
window (76 classes, commuting probability 0.006); the non-metacyclic enumeration by block decompositions is in
`../lift_nonmetacyclic_groups/`.
