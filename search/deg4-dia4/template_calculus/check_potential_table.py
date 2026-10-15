"""Independent re-check of the potential table (potential3.py): for every connector type
x in Z^3 with |x| <= 3 (degree 4 - |x| >= 1; 20 types) and every unordered pair, E(x, x') =
2 v + node(x) + node(x') <= 0, with the zero set exactly {complementary triple pairs} u
{d_ij - s_k, {i,j,k} = {1,2,3}}; quads (|x| = 4, isolated) have self - 15 < 0. Exact rationals."""
from fractions import Fraction as Fr
import itertools, json

types = [x for x in itertools.product(range(5), repeat=3) if sum(x) <= 3]
assert len(types) == 20
def self_(x): return sum(x) ** 2 - sum(t * t for t in x)
def node(x): return Fr(self_(x) - 3 * sum(x) - 3, 4 - sum(x))
def v(x, y): return sum(x) * sum(y) - sum(a * b for a, b in zip(x, y))
zero, pos = [], []
for x, y in itertools.combinations_with_replacement(types, 2):
    E = 2 * v(x, y) + node(x) + node(y)
    if E > 0: pos.append((x, y, str(E)))
    if E == 0: zero.append((x, y))
def is_unit(x, y): return sum(x) == 3 and sum(y) == 3 and tuple(a + b for a, b in zip(x, y)) == (2, 2, 2)
def is_ds(x, y):
    for d, s in ((x, y), (y, x)):
        if sum(d) == 2 and max(d) == 1 and sum(s) == 2 and max(s) == 2 and d[s.index(2)] == 0:
            return True
    return False
expected = [(x, y) for x, y in itertools.combinations_with_replacement(types, 2) if is_unit(x, y) or is_ds(x, y)]
quads = {x: self_(x) - 15 for x in itertools.product(range(5), repeat=3) if sum(x) == 4}
out = {"types": len(types), "pairs": 210, "positive_E": pos, "zero_E": zero, "zero_E_as_expected": sorted(zero) == sorted(expected),
       "n_zero": len(zero), "quad_potentials": sorted(set(quads.values())), "all_quads_negative": all(p < 0 for p in quads.values())}
print(json.dumps(out, indent=1, default=str))
