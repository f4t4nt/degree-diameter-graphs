"""Convert the text output of clusters108.c into the jsonl format of clusters108.py
(fields W, deficit, cost8, slots, G, quads, ex, hh, nres, decl, sv22) plus the .meta.json window
table, so that pairs108.py / combine108.py / windows108.py shapes can consume it unchanged.

Usage: c2jsonl.py CMAX IN.txt OUT.jsonl
"""
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import clusters108 as K  # noqa: E402


def main(args):
    cmax = int(args[0])
    inp, out = args[1], args[2]
    wins = K.window_table(cmax * 8, False, sys.stderr)
    n = 0
    with open(out, "w") as fout:
        for line in open(inp):
            p = line.split()
            if not p or p[0] != "W":
                continue
            ws = []
            i = 1
            while p[i] != "D":
                a, b = p[i].split(":")
                ws.append((int(a), int(b)))
                i += 1
            decl, nres = int(p[i + 1]), int(p[i + 3])
            a = [[0] * K.NT for _ in range(K.NT)]
            sv = [K.collections.Counter() for _ in range(3)]
            slots, ex, hh = [0, 0, 0], [0] * 6, [0, 0, 0]
            st = {"deficit": 0, "cost8": 0, "G": 0, "quads": 0}
            for wid, cnt in ws:
                K.Growth.add(wins[wid], cnt, a, sv, slots, ex, hh, st)
            rec = {"W": ws, "deficit": st["deficit"], "cost8": st["cost8"], "slots": slots, "G": st["G"], "quads": st["quads"],
                   "ex": ex, "hh": hh, "nres": nres, "decl": decl, "sv22": [sv[j].get((2, 2), 0) for j in range(3)]}
            fout.write(json.dumps(rec) + "\n")
            n += 1
    summary = {"source": inp, "clusters": n, "cmax": cmax, "windows": len(wins), "note": "clusters108.c output converted by c2jsonl.py"}
    json.dump({"summary": summary, "windows": [[list(w.x), [list(t) for t in w.M]] for w in wins]}, open(out + ".meta.json", "w"))
    print(f"{n} clusters -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1:])
