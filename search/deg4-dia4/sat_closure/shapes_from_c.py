"""Finish the combine of combine108.c in Python: read its survivors ("u G cost8 nparts | wid:cnt ...
| parts ..."), deduplicate the window multisets under block relabelling (windows108.canon_W), apply
the slot-vector matching relaxation (shapes107.matching_ok on the one-connector-per-slot gadget)
and write the window multisets in the jsonl format of pairs108.py / windows108.py shapes.

Usage: shapes_from_c.py CMAX IN.txt OUT.jsonl
"""
import collections
import hashlib
import json
import pathlib
import sys
from fractions import Fraction as FR

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import clusters108 as K  # noqa: E402
import windows108 as W  # noqa: E402


def main(args):
    cmax = int(args[0])
    inp, out = args[1], args[2]
    wins = K.window_table(cmax * 8, False, sys.stderr)
    wid_of = {(w.x, w.M): w.id for w in wins}
    stats = collections.Counter()
    seen = set()
    n = 0
    with open(out, "w") as fout:
        for line in open(inp):
            head, body, parts = line.split("|")
            u, G, cost8, nparts = (int(v) for v in head.split())
            items = sorted((int(a), int(b)) for a, b in (p.split(":") for p in body.split()))
            stats["records"] += 1
            key = W.canon_W(items, wins, wid_of)
            dig = hashlib.blake2b(repr(key).encode(), digest_size=16).digest()
            if dig in seen:
                stats["relabel_dup"] += 1
                continue
            seen.add(dig)
            svl = []
            for wid, cnt in items:
                svl.extend(wins[wid].svec * cnt)
            if not W.S.matching_ok([W.sv_gadget(svl)]):
                stats["sv_matching_fail"] += 1
                continue
            stats["survive"] += 1
            stats[f"survive_cost_{cost8 // 8}"] += 1
            stats[f"survive_parts_{nparts}"] += 1
            rec = {"u": u, "G": G, "cost": cost8 // 8, "phi": str(-FR(cost8, 8)), "clusters": nparts, "parts": parts.split()[1:],
                   "W": [[list(wins[wid].x), [list(t) for t in wins[wid].M], cnt] for wid, cnt in items]}
            fout.write(json.dumps(rec) + "\n")
            n += 1
    print(json.dumps({"in": inp, "out": out, "stats": dict(stats)}), file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1:])
