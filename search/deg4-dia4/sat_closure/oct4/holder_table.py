"""Holder table for the Oct x4 tight structure (RESULTS.md section 3), computed from the component
definitions of the structure theorem rather than by hand.

Component kinds (potential.py): unit = two triples with profile sum s, a permutation of (2,2,1,1);
cycle double of type {j,k} whose two neighbours have the complementary type {l,m}; same-block
double 2 e_m whose neighbours are straight triples missing block m; straight triple missing m
adjacent to such a double. For every kind and every block j in which the holder has a port, the
script lists the carried vector f restricted to the other blocks and its total g, then all pairs
of holders of a common port of block j whose vectors add to (3,3,3) as the holder lemma requires.
"""
import itertools
import json
import sys

NB = 4


def e(*blocks):
    return tuple(1 if b in blocks else 0 for b in range(NB))


def add(*vs):
    return tuple(sum(t) for t in zip(*vs))


def holders():
    """(description, own profile x, carried vector f over all blocks)."""
    out = []
    for heavy in itertools.combinations(range(NB), 2):
        s = tuple(2 if b in heavy else 1 for b in range(NB))
        for x in itertools.product(range(3), repeat=NB):
            if sum(x) != 3 or any(x[b] > s[b] for b in range(NB)):
                continue
            out.append((f"unit{heavy} triple x={x}", x, s))
    for t in itertools.combinations(range(NB), 2):
        comp = tuple(b for b in range(NB) if b not in t)
        out.append((f"cycle double {t}", e(*t), add(e(*t), e(*comp), e(*comp))))
    for m in range(NB):
        others = tuple(b for b in range(NB) if b != m)
        s = add(e(m), e(m))
        out.append((f"same-block double block {m}", s, add(s, e(*others), e(*others))))
        out.append((f"straight triple missing {m} next to same-block double", e(*others), add(e(*others), s)))
    return out


def main():
    H = holders()
    rows = []
    for name, x, f in H:
        for j in range(NB):
            if x[j] == 0:
                continue
            foreign = tuple(f[i] for i in range(NB) if i != j)
            rows.append({"holder": name, "port_block": j, "f_foreign": foreign, "g": sum(foreign)})
    gs = sorted(set(r["g"] for r in rows))
    pairs = []
    for r1, r2 in itertools.combinations_with_replacement(rows, 2):
        if r1["port_block"] != r2["port_block"]:
            continue
        if all(a + b == 3 for a, b in zip(r1["f_foreign"], r2["f_foreign"])):
            pairs.append((r1["holder"], r2["holder"], r1["port_block"]))
    def kind(name):
        if name.startswith("unit"):
            return name.split(" triple")[0]
        if name.startswith("straight"):
            return "path triple (next to same-block double)"
        if name.startswith("same-block"):
            return "same-block double"
        return name

    kinds = sorted(set((kind(p[0]), kind(p[1])) for p in pairs))
    rec = {"holder_rows": len(rows), "g_values": gs,
           "g_by_kind": sorted(set((r["holder"].split(" x=")[0] if "unit" not in r["holder"] else
                                    ("unit heavy" if r["f_foreign"].count(2) == 1 else "unit light"), r["g"]) for r in rows)),
           "valid_pairs": len(pairs), "valid_pair_kinds": kinds,
           "same_block_double_in_valid_pair": any(kind(p[0]) == "same-block double" or kind(p[1]) == "same-block double" for p in pairs),
           "two_doubles_in_valid_pair": any(p[0].startswith("cycle") and p[1].startswith("cycle") for p in pairs),
           "pairs": pairs}
    print(json.dumps({k: rec[k] for k in rec if k != "pairs"}, indent=1), file=sys.stderr)
    json.dump(rec, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
