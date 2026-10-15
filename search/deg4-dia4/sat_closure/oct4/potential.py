"""Potential-function proof of the Oct x4 structure theorem (finite table check).

Notation: connector c has port profile x = x(c) in Z^4 (ports per block), connector degree d(c),
|x| + d <= 4. T(c) = |x||F| - x.F with F = x + sum of neighbours' x, so that the hub-to-port
identity reads 2 sum_c T(c) >= 864, i.e. sum_c T(c) >= 432, and
    T(c) = self(c) + sum_{c' ~ c} v(c, c'),  self(x) = |x|^2 - |x|_2^2,  v(x, x') = |x||x'| - x.x'.

Potential: phi(c) = T(c) - 3|x| - 4. With 36 connectors and at most 96 filled slots,
    sum_c phi(c) = sum_c T(c) - 3 * slots - 144 >= 432 - 288 - 144 = 0.
Distribute phi over edges: for an edge cc' let
    E(c, c') = 2 v(c, c') + (self(x) - 3|x| - 4)/d(c) + (self(x') - 3|x'| - 4)/d(c'),
then sum_c phi(c) = sum_{edges} E(c, c') + sum_{isolated c} (self(x) - 3|x| - 4).
This script checks over all types (x, d) with d >= 1 and all pairs that E <= 0, lists the pairs with
E = 0, and checks that isolated connectors (d = 0) have negative potential. Consequently
sum_c phi(c) <= 0, so equality holds throughout: all 96 slots are filled, no connector is
isolated, and every edge is one of the E = 0 types. Expected E = 0 types: two triples whose
profiles add to a permutation of (2,2,1,1); two cross doubles on complementary block pairs; a
same-block double and a straight triple missing its block.
"""
import itertools
import json
import sys
from fractions import Fraction

NB = 4


def types():
    out = []
    for x in itertools.product(range(5), repeat=NB):
        if sum(x) > 4:
            continue
        for d in range(0, 5 - sum(x)):
            out.append((x, d))
    return out


def selfv(x):
    return sum(x) ** 2 - sum(a * a for a in x)


def v(x, u):
    return sum(x) * sum(u) - sum(a * b for a, b in zip(x, u))


def node_term(x, d):
    return Fraction(selfv(x) - 3 * sum(x) - 4, d)


def main():
    TY = types()
    iso = [(t, selfv(t[0]) - 3 * sum(t[0]) - 4) for t in TY if t[1] == 0]
    assert all(val < 0 for _, val in iso), iso
    FREE = [t for t in TY if t[1] >= 1]
    zero, worst = [], Fraction(-100)
    npairs = 0
    for i, t in enumerate(FREE):
        for u in FREE[i:]:
            e = 2 * v(t[0], u[0]) + node_term(*t) + node_term(*u)
            npairs += 1
            assert e <= 0, (t, u, e)
            if e == 0:
                zero.append((t, u))
            elif e > worst:
                worst = e
    # classify the zero-E pairs
    classes = {"unit": 0, "complementary_cross_doubles": 0, "same_block_double_straight_triple": 0, "other": []}
    for t, u in zero:
        x, dx = t
        y, dy = u
        s = tuple(a + b for a, b in zip(x, y))
        if sum(x) == 3 and sum(y) == 3 and dx == dy == 1 and sorted(s) == [1, 1, 2, 2]:
            classes["unit"] += 1
        elif sum(x) == 2 and sum(y) == 2 and dx == dy == 2 and max(x) == 1 and max(y) == 1 and v(x, y) == 4:
            classes["complementary_cross_doubles"] += 1
        elif dx == dy and False:
            pass
        elif {sum(x), sum(y)} == {2, 3}:
            dbl, tri = (t, u) if sum(x) == 2 else (u, t)
            if dbl[1] == 2 and tri[1] == 1 and max(dbl[0]) == 2 and max(tri[0]) == 1 and v(dbl[0], tri[0]) == 6:
                classes["same_block_double_straight_triple"] += 1
            else:
                classes["other"].append((t, u))
        else:
            classes["other"].append((t, u))
    rec = {"types": len(TY), "isolated_types": len(iso), "edge_type_pairs": npairs, "zero_E_pairs": len(zero),
           "largest_negative_E": str(worst), "classes": {k: (val if k != "other" else val) for k, val in classes.items()},
           "zero_E_list": [[list(t[0]), t[1], list(u[0]), u[1]] for t, u in zero]}
    print(json.dumps({k: rec[k] for k in rec if k != "zero_E_list"}), file=sys.stderr)
    json.dump(rec, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
