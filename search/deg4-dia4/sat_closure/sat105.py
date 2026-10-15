"""Exact completion search over connector splits for one tight unit design (skeleton).

Input: a skeleton = 16 units, each owning two ports per block (runs/unit_design_classes.json).
Every unit is realised either as a matched pair of triples (10 splits: connector x holds the unit's
first block-1 port plus two of the other five ports, y holds the rest) or, for exactly one unit in
the 105 case, as a triangle gadget of three doubles:
  (a) cross triangle   x = d12 (p, q), y = d23 (q', r), z = d13 (p', r')   8 orientations
  (b) same-block       x = s1 (p, p'), y = s2 (q, q'), z = s3 (r, r')      1
  (c) two cross + one  x = d_ij (p, q), y = d_ij' (p', q'), z = s_k (r, r')  6 (block k and pairing)
The 104 mode allows no gadget and tests which skeletons admit a diameter-4 realisation at all.

Encoding (CaDiCaL through python-sat): option selectors per unit with exactly-one; adjacency
literals defined from selectors; for every vertex pair not automatically within distance 4 (all
within-block pairs are), a clause over witness walks of length <= 4 in the potential graph, each
witness a conjunction of the variable adjacency literals on the walk. Solutions are decoded and
verified independently by BFS before being reported.
"""
import collections
import hashlib
import itertools
import json
import sys
import time

from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Cadical153


import os

# the 54 skeletons of the template calculus (output of template_calculus/unit_designs2.py)
SKELETON_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "skeletons54.json")


def load_skeletons():
    data = json.load(open(SKELETON_FILE))
    skels = []
    for c in data["classes"]:
        units = []
        for t in c["design"]:
            per_block = []
            for b in range(3):
                edges = EDGES[t[b]]
                per_block.append([("P", b, l, r) for (l, r) in edges])
            units.append(per_block)
        skels.append(units)
    return skels


EDGES = {
    "A1": ((0, 0), (1, 1)), "A2": ((0, 1), (1, 0)),
    "B1": ((2, 2), (3, 3)), "B2": ((2, 3), (3, 2)),
    "C1": ((0, 2), (1, 3)), "C2": ((0, 3), (1, 2)),
    "D1": ((2, 0), (3, 1)), "D2": ((2, 1), (3, 0)),
}


def unit_options(ports, allow_gadget):
    """ports: [[p, p'], [q, q'], [r, r']]. Returns list of options:
    (kind, holds: {0: set, 1: set, 2: set}, extra_edges: set of (i, j) among {0,1,2}, active_z)."""
    (p, p2), (q, q2), (r, r2) = ports
    opts = []
    others = [p2, q, q2, r, r2]
    for extra in itertools.combinations(others, 2):
        hx = {p} | set(extra)
        hy = set([p, p2, q, q2, r, r2]) - hx
        opts.append(("pair", {0: hx, 1: hy, 2: set()}, set(), False))
    if allow_gadget:
        for sp, sq, sr in itertools.product((0, 1), repeat=3):
            pp, pp2 = (p, p2) if sp == 0 else (p2, p)
            qq, qq2 = (q, q2) if sq == 0 else (q2, q)
            rr, rr2 = (r, r2) if sr == 0 else (r2, r)
            opts.append(("tri-a", {0: {pp, qq}, 1: {qq2, rr}, 2: {pp2, rr2}}, {(1, 2), (0, 2)}, True))
        opts.append(("tri-b", {0: {p, p2}, 1: {q, q2}, 2: {r, r2}}, {(1, 2), (0, 2)}, True))
        for k, (X, Y, Zs) in enumerate((((p, p2), (q, q2), (r, r2)), ((q, q2), (r, r2), (p, p2)), ((p, p2), (r, r2), (q, q2)))):
            for sw in (0, 1):
                a, a2 = X
                b, b2 = Y if sw == 0 else (Y[1], Y[0])
                opts.append(("tri-c", {0: {a, b}, 1: {a2, b2}, 2: set(Zs)}, {(1, 2), (0, 2)}, True))
    return opts


class Encoder:
    def __init__(self, units, mode):
        self.units = units
        self.mode = mode  # "104" or "105"
        self.pool = IDPool()
        self.cnf = CNF()
        self.hubs = [("H", b, side, i) for b in range(3) for side in ("L", "R") for i in range(4)]
        self.ports = [("P", b, l, r) for b in range(3) for l in range(4) for r in range(4)]
        self.conns = [("C", u, k) for u in range(len(units)) for k in range(3)]
        self.block = {}
        for v in self.hubs + self.ports:
            self.block[v] = v[1]
        self.fixed = collections.defaultdict(set)   # fixed adjacency
        for (_, b, l, r) in self.ports:
            self.fixed[("P", b, l, r)] |= {("H", b, "L", l), ("H", b, "R", r)}
            self.fixed[("H", b, "L", l)].add(("P", b, l, r))
            self.fixed[("H", b, "R", r)].add(("P", b, l, r))
        self.options = {}
        self.sel = {}
        self.adj_lit = {}       # (connector, port) -> var ; (connector, connector) -> var or True
        self.active = {}
        self.build_units()

    def var(self, name):
        return self.pool.id(name)

    def build_units(self):
        gadget_lits = []
        for u, ports in enumerate(self.units):
            opts = unit_options(ports, self.mode in ("105", "120"))
            self.options[u] = opts
            lits = [self.var(("sel", u, o)) for o in range(len(opts))]
            self.sel[u] = lits
            for cl in CardEnc.equals(lits=lits, bound=1, vpool=self.pool, encoding=EncType.pairwise).clauses:
                self.cnf.append(cl)
            # activity of z
            za = self.var(("active", u))
            self.active[("C", u, 2)] = za
            gl = [lits[o] for o, opt in enumerate(opts) if opt[3]]
            # za <-> OR(gl)
            if gl:
                self.cnf.append([-za] + gl)
                for g in gl:
                    self.cnf.append([za, -g])
                gadget_lits += gl
            else:
                self.cnf.append([-za])
            allports = [t for blk in ports for t in blk]
            for k in range(3):
                cv = ("C", u, k)
                for t in allports:
                    holders = [lits[o] for o, opt in enumerate(opts) if t in opt[1][k]]
                    a = self.var(("adj", cv, t))
                    self.adj_lit[(cv, t)] = a
                    self.adj_lit[(t, cv)] = a
                    if holders:
                        self.cnf.append([-a] + holders)
                        for h in holders:
                            self.cnf.append([a, -h])
                    else:
                        self.cnf.append([-a])
            # connector-connector edges
            self.adj_lit[(("C", u, 0), ("C", u, 1))] = True
            self.adj_lit[(("C", u, 1), ("C", u, 0))] = True
            for (i, j) in ((1, 2), (0, 2)):
                e = self.var(("cedge", u, i, j))
                self.adj_lit[(("C", u, i), ("C", u, j))] = e
                self.adj_lit[(("C", u, j), ("C", u, i))] = e
                gl2 = [lits[o] for o, opt in enumerate(opts) if (i, j) in opt[2]]
                if gl2:
                    self.cnf.append([-e] + gl2)
                    for g in gl2:
                        self.cnf.append([e, -g])
                else:
                    self.cnf.append([-e])
        # gadget count: 105 = exactly one unit is a gadget; 120 = every unit is a gadget
        if self.mode == "105":
            for cl in CardEnc.equals(lits=gadget_lits, bound=1, vpool=self.pool, encoding=EncType.seqcounter).clauses:
                self.cnf.append(cl)
        elif self.mode == "120":
            for u in range(len(self.units)):
                for o, opt in enumerate(self.options[u]):
                    if not opt[3]:
                        self.cnf.append([-self.sel[u][o]])
        # potential adjacency
        self.pot = collections.defaultdict(set)
        for v, s in self.fixed.items():
            self.pot[v] |= s
        for (a, b), lit in self.adj_lit.items():
            self.pot[a].add(b)

    def lit_of_edge(self, a, b):
        if b in self.fixed.get(a, ()):
            return True
        return self.adj_lit.get((a, b))

    def witnesses(self, src, maxlen=4):
        """For each vertex w reachable within maxlen from src in the potential graph, the set of
        minimal literal-sets of walks (True = unconditional)."""
        found = collections.defaultdict(set)   # w -> set of frozenset(lits)
        uncond = set()
        # BFS over (vertex, literal set) with pruning of supersets is expensive; enumerate walks
        stack = [(src, frozenset(), 0, None)]
        while stack:
            v, lits, d, prev = stack.pop()
            if v != src:
                if not lits:
                    uncond.add(v)
                elif v not in uncond:
                    found[v].add(lits)
            if d == maxlen:
                continue
            for w in self.pot[v]:
                if w == prev:
                    continue
                l = self.lit_of_edge(v, w)
                if l is None:
                    continue
                nl = lits if l is True else lits | {l}
                if len(nl) > 4:
                    continue
                stack.append((w, nl, d + 1, v))
        for v in uncond:
            found.pop(v, None)
        return uncond, found

    def encode_distances(self):
        verts = self.hubs + self.ports + self.conns
        index = {v: i for i, v in enumerate(verts)}
        n_pairs = 0
        n_wit = 0
        auto = 0
        for i, v in enumerate(verts):
            uncond, found = self.witnesses(v)
            for w in verts[i + 1:]:
                if v in self.block and w in self.block and self.block[v] == self.block[w]:
                    continue  # within-block pairs are within 4 by block edges
                if w in uncond:
                    auto += 1
                    continue
                sets = found.get(w, set())
                # minimalise
                sets = sorted(sets, key=len)
                minimal = []
                for s in sets:
                    if not any(m <= s for m in minimal):
                        minimal.append(s)
                clause = []
                cond = []
                for ep in (v, w):
                    if ep[0] == "C" and ep[2] == 2:
                        cond.append(-self.active[ep])
                for s in minimal:
                    if len(s) == 1:
                        clause.append(next(iter(s)))
                    else:
                        wv = self.var(("wit", v, w, tuple(sorted(s))))
                        for l in s:
                            self.cnf.append([-wv, l])
                        clause.append(wv)
                    n_wit += 1
                self.cnf.append(cond + clause)
                n_pairs += 1
        return n_pairs, n_wit, auto

    def decode(self, model):
        pos = set(l for l in model if l > 0)
        edges = set()
        for v, s in self.fixed.items():
            for w in s:
                edges.add(tuple(sorted((str(v), str(w)))))
        active = set()
        for (a, b), lit in self.adj_lit.items():
            if lit is True or lit in pos:
                edges.add(tuple(sorted((str(a), str(b)))))
        for u in range(len(self.units)):
            if self.active[("C", u, 2)] in pos:
                active.add(u)
        chosen = {u: next(o for o, l in enumerate(self.sel[u]) if l in pos) for u in range(len(self.units))}
        return sorted(edges), active, chosen


def verify_edges(edges):
    adj = collections.defaultdict(set)
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    verts = sorted(adj)
    degs = collections.Counter(len(adj[v]) for v in verts)
    bad = 0
    diam = 0
    for v in verts:
        dist = {v: 0}
        fr = [v]
        while fr:
            nx_ = []
            for x in fr:
                for y in adj[x]:
                    if y not in dist:
                        dist[y] = dist[x] + 1
                        nx_.append(y)
            fr = nx_
        if len(dist) != len(verts):
            bad += len(verts) - len(dist)
        diam = max(diam, max(dist.values()))
        bad += sum(1 for u, d in dist.items() if d > 4 and u > v)
    return {"n": len(verts), "degrees": dict(degs), "diameter": diam, "bad_pairs": bad}


def run(skel_index, units, mode, time_limit):
    t0 = time.time()
    enc = Encoder(units, mode)
    n_pairs, n_wit, auto = enc.encode_distances()
    build = time.time() - t0
    rec = {"skeleton": skel_index, "mode": mode, "vars": enc.pool.top, "clauses": len(enc.cnf.clauses),
           "pairs_encoded": n_pairs, "witnesses": n_wit, "auto_pairs": auto, "build_s": round(build, 1)}
    solver = Cadical153(bootstrap_with=enc.cnf.clauses)
    t1 = time.time()
    sat = solver.solve()
    rec["solve_s"] = round(time.time() - t1, 1)
    rec["sat"] = sat
    if sat:
        edges, active, chosen = enc.decode(solver.get_model())
        v = verify_edges(edges)
        rec["verified"] = v
        rec["gadget_units"] = sorted(active)
        rec["chosen_option_kinds"] = collections.Counter(enc.options[u][o][0] for u, o in chosen.items())
        rec["edges_sha256"] = hashlib.sha256(json.dumps(edges).encode()).hexdigest()
        with open(f"runs/sat105_solution_skel{skel_index}_{mode}.json", "w") as f:
            json.dump({"edges": edges, "gadget_units": sorted(active), "verified": v}, f)
    solver.delete()
    return rec


def main():
    mode = sys.argv[1]
    which = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 and sys.argv[2] != "all" else None
    skels = load_skeletons()
    out = []
    for i, units in enumerate(skels):
        if which is not None and i not in which:
            continue
        rec = run(i, units, mode, None)
        print(json.dumps(rec, default=str), file=sys.stderr, flush=True)
        out.append(rec)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
