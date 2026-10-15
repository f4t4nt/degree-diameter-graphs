"""Export the window table (windows108.py / clusters108.py, committed-cost filter)
as a plain text file for clusters108.c, so that the C enumeration uses exactly the same windows,
ids, costs and derived tables as the validated Python growth.

Usage: export_windows.py CMAX --out=runs/windows_cN.txt     (CMAX = max cost, e.g. 12 or 6)

Format (all integers; costs in eighths):
  NT 20
  TYPES  (20 lines: t x0 x1 x2 deg def)
  HC     (20 lines of 20 ints: HC[t][m] = cost in eighths of one half-edge (m -> t))
  NW n
  n window lines: id t nx quad x0 x1 x2 cost8 deficit nh [m k]*nh ns [j f0 f1]*ns ex0..ex5 hh0 hh1 hh2
  (t = -1 for quads; halves (m, k) = neighbour type m with multiplicity k; svec entries with
  multiplicity, one line per window in id order.)
"""
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import clusters108 as K  # noqa: E402
import windows108 as W  # noqa: E402


def main(args):
    cmax = int(args[0])
    out = next(a.split("=")[1] for a in args if a.startswith("--out="))
    wins = K.window_table(cmax * 8, False, sys.stderr)
    lines = [f"NT {K.NT}"]
    for t, x in enumerate(W.NONISO):
        lines.append(f"{t} {x[0]} {x[1]} {x[2]} {K.DEG[t]} {K.DEF[t]}")
    for t in range(K.NT):
        lines.append(" ".join(str(K.HC[t][m]) for m in range(K.NT)))
    lines.append(f"NW {len(wins)}")
    for w in wins:
        parts = [w.id, -1 if w.quad else w.t, w.nx, int(w.quad), *w.x, w.cost8, w.deficit, len(w.halves)]
        for m, k in w.halves:
            parts += [m, k]
        parts.append(len(w.svec))
        for (j, fv) in w.svec:
            parts += [j, fv[0], fv[1]]
        parts += list(w.ex) + list(w.hh)
        lines.append(" ".join(str(v) for v in parts))
    pathlib.Path(out).write_text("\n".join(lines) + "\n")
    print(f"{len(wins)} windows written to {out}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1:])
