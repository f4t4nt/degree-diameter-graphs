"""Write the 36 diameter-4 graphs on 104 vertices as edge lists in the block numbering.

Input: the class files written by sat_closure/classify104.py (one JSON list per skeleton, each
entry {"edges": [[u, v], ...], "known_index": k or null} with the symbolic labels of the encoder:
("H", block, side, i) hubs, ("P", block, l, r) ports, ("C", unit, member) connectors).

Output: classNN.edges.txt with the integer numbering

    block b in {0, 1, 2} occupies 24 b .. 24 b + 23:
        hubs L_i = 24 b + i, R_i = 24 b + 4 + i (i = 0..3),
        port P_xy between L_x and R_y = 24 b + 8 + 4 x + y;
    connector pair k (k = 0..15) is 72 + 2 k and 73 + 2 k.

NN is the known_index of the entry (the order of the campaign's corpus of previously found
graphs) and the next free number for entries without one.

Usage: python export_classes.py OUTDIR CLASSES.json [CLASSES.json ...]
"""
import ast
import json
import pathlib
import sys


def label_to_int(label):
    t = ast.literal_eval(label)
    if t[0] == "H":
        _, b, side, i = t
        return 24 * b + (0 if side == "L" else 4) + i
    if t[0] == "P":
        _, b, l, r = t
        return 24 * b + 8 + 4 * l + r
    if t[0] == "C":
        _, k, m = t
        return 72 + 2 * k + m
    raise ValueError(label)


def kind(v):
    if v >= 72:
        return "connector"
    return "hub" if v % 24 < 8 else "port"


def write_class(path, nn, edges):
    by_type = {"hub-port edges, block A": [], "hub-port edges, block B": [], "hub-port edges, block C": [],
               "port-connector edges": [], "connector-connector edges": []}
    for u, v in edges:
        ku, kv = kind(u), kind(v)
        if {ku, kv} == {"hub", "port"}:
            by_type["hub-port edges, block " + "ABC"[u // 24]].append((u, v))
        elif {ku, kv} == {"port", "connector"}:
            by_type["port-connector edges"].append((u, v))
        elif ku == kv == "connector":
            by_type["connector-connector edges"].append((u, v))
        else:
            raise ValueError((u, v, ku, kv))
    lines = [f"# deg 4 diam 4 n 104 class {nn:02d}",
             "# blocks A = 0..23, B = 24..47, C = 48..71: hubs L0-L3 = +0..3, R0-R3 = +4..7, port P_xy = +8 + 4x + y",
             "# connectors 72..103 in matched pairs (72 + 2k, 73 + 2k), k = 0..15"]
    for section, es in by_type.items():
        lines.append(f"# {section} ({len(es)})")
        lines.extend(f"{u} {v}" for u, v in sorted(es))
    path.write_text("\n".join(lines) + "\n")


def main():
    outdir = pathlib.Path(sys.argv[1])
    outdir.mkdir(parents=True, exist_ok=True)
    entries = []
    for f in sys.argv[2:]:
        entries.extend(json.load(open(f)))
    known = {e["known_index"]: e for e in entries if e.get("known_index") is not None}
    unknown = [e for e in entries if e.get("known_index") is None]
    numbered = [(k, known[k]) for k in sorted(known)]
    nxt = max(known) + 1 if known else 0
    for e in unknown:
        numbered.append((nxt, e))
        nxt += 1
    for nn, e in numbered:
        edges = sorted(tuple(sorted((label_to_int(u), label_to_int(v)))) for u, v in e["edges"])
        assert len(edges) == 208 and len(set(edges)) == 208
        write_class(outdir / f"class{nn:02d}.edges.txt", nn, edges)
    print(f"wrote {len(numbered)} classes to {outdir}")


if __name__ == "__main__":
    main()
