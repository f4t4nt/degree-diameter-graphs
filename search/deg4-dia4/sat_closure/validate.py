"""Validation of the SAT completion encoder (sat105.Encoder).

1. Positive controls: for each of the known 104-vertex classes (../classes/), build the encoder on
   that graph's own unit design (its real blocks and ports, hubs relabelled 0..3 per side), fix every
   unit's option to the split the known graph actually uses, and check satisfiability. A sound encoder
   must say SAT for all of them: these graphs have diameter 4, so no clause may exclude them.
2. Gadget positive control: on one skeleton in 105 mode, drop every distance clause that involves a
   gadget vertex (connectors z, and x/y of the gadget unit are unknown in advance, so drop all clauses
   mentioning connector vertices of units when the unit is the gadget, by encoding only non-connector
   pairs) and check the solver can place a triangle; decode and check the gadget wiring and degrees.
   This exercises the gadget option encoding on a satisfiable instance.
"""
import collections
import json
import pathlib
import sys

from pysat.solvers import Cadical153

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "template_calculus"))
import sat105  # noqa: E402
import surgery  # noqa: E402


def known_encoding_check(name, verts, adj, trip, blocks):
    T = surgery.Template(verts, adj, trip, blocks)
    # hub relabelling per block: L hubs and R hubs to 0..3
    hub_idx = {}
    for b, blk in enumerate(trip):
        L, R, pmap = blocks[blk]
        for i, h in enumerate(L):
            hub_idx[h] = ("L", i)
        for i, h in enumerate(R):
            hub_idx[h] = ("R", i)
    def port_key(q):
        l, r = T.port_hubs[q]
        b = T.port_block[q]
        li = hub_idx[l][1] if hub_idx[l][0] == "L" else hub_idx[r][1]
        ri = hub_idx[r][1] if hub_idx[r][0] == "R" else hub_idx[l][1]
        return ("P", b, li, ri)
    units = []
    unit_pairs = T.pairs
    for pair in unit_pairs:
        per = T.pair_ports[pair]
        units.append([[port_key(q) for q in per[b]] for b in range(3)])
    enc = sat105.Encoder(units, "104")
    enc.encode_distances()
    assumptions = []
    for u, pair in enumerate(unit_pairs):
        c, c2 = pair
        ports_c = {port_key(q) for q in T.ports_of[c]}
        ports_c2 = {port_key(q) for q in T.ports_of[c2]}
        first = units[u][0][0]
        xset = ports_c if first in ports_c else ports_c2
        opts = enc.options[u]
        idx = [o for o, opt in enumerate(opts) if opt[1][0] == xset]
        assert len(idx) == 1, (u, xset, [opt[1][0] for opt in opts])
        assumptions.append(enc.sel[u][idx[0]])
    solver = Cadical153(bootstrap_with=enc.cnf.clauses)
    sat = solver.solve(assumptions=assumptions)
    rec = {"class": name, "fixed_split_sat": sat}
    if sat:
        edges, active, chosen = enc.decode(solver.get_model())
        rec["verified"] = sat105.verify_edges(edges)
    solver.delete()
    return rec


def gadget_control():
    units = sat105.load_skeletons()[50]
    enc = sat105.Encoder(units, "105")
    # encode only pairs among hubs and ports (no connector endpoints) so the instance is satisfiable
    verts = enc.hubs + enc.ports
    n_pairs = 0
    for i, v in enumerate(verts):
        uncond, found = enc.witnesses(v)
        for w in verts[i + 1:]:
            if enc.block[v] == enc.block[w] or w in uncond:
                continue
            sets = sorted(found.get(w, set()), key=len)
            minimal = []
            for s in sets:
                if not any(m <= s for m in minimal):
                    minimal.append(s)
            clause = []
            for s in minimal:
                if len(s) == 1:
                    clause.append(next(iter(s)))
                else:
                    wv = enc.var(("wit", v, w, tuple(sorted(s))))
                    for l in s:
                        enc.cnf.append([-wv, l])
                    clause.append(wv)
            enc.cnf.append(clause)
            n_pairs += 1
    solver = Cadical153(bootstrap_with=enc.cnf.clauses)
    sat = solver.solve()
    rec = {"gadget_control_pairs_encoded": n_pairs, "sat": sat}
    if sat:
        edges, active, chosen = enc.decode(solver.get_model())
        v = sat105.verify_edges(edges)
        rec["verified"] = v
        rec["gadget_units"] = sorted(active)
        rec["gadget_kind"] = [enc.options[u][chosen[u]][0] for u in active]
        # hub/port pairs must all be within 4; connector pairs may fail (not encoded)
        adjd = collections.defaultdict(set)
        for a, b in edges:
            adjd[a].add(b)
            adjd[b].add(a)
        hubports = [str(x) for x in verts]
        bad_hp = 0
        for s in hubports:
            dist = {s: 0}
            fr = [s]
            while fr:
                nx_ = []
                for x in fr:
                    for y in adjd[x]:
                        if y not in dist:
                            dist[y] = dist[x] + 1
                            nx_.append(y)
                fr = nx_
            bad_hp += sum(1 for t in hubports if dist.get(t, 99) > 4)
        rec["hub_port_pairs_beyond_4"] = bad_hp // 2
    solver.delete()
    return rec


def main():
    out = {"known": [], "gadget_control": None}
    for cls in surgery.load_classes():
        rec = known_encoding_check(*cls)
        print(rec, file=sys.stderr, flush=True)
        out["known"].append(rec)
    out["known_all_sat"] = all(r["fixed_split_sat"] for r in out["known"])
    out["gadget_control"] = gadget_control()
    print(out["gadget_control"], file=sys.stderr, flush=True)
    json.dump(out, sys.stdout, indent=1, default=str)


if __name__ == "__main__":
    main()
