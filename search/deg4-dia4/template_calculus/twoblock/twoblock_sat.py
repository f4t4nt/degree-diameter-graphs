"""SAT model for the two-block architecture with a FREE remainder W.

Two vertex-disjoint subdivided K4,4 blocks (hubs saturated, no port-port edge) are fixed; the
remainder W has m = n - 48 vertices, and both its plain edges E[u][v] and its port attachments
A[c][q] are free variables (4-regularity: every W-vertex has degree 4 counting ports, every port has
two holders). Diameter 4 is encoded for every pair not automatically within 4 through
N2[u][v] ("walk of length <= 2", defined one-directionally by its witnesses) and the identity
d(u, v) <= 4  <=>  exists w: N2[u][w] and N2[w][v].

Families (restrictions on W, all optional):
  general  no restriction (any profiles, any plain graph)
  F2       W[0:32] hold exactly two ports each (doubles), W[32:] hold none (routers)
  F2M      F2 and the doubles form 16 adjacent partner pairs (2k, 2k+1), each double having
           exactly one further plain edge (to a router): the record's connector architecture with
           the third block replaced by an arbitrary router graph; the facts forced by B = 256
           (RESULTS.md section 2: unit profile (2,2), 2-matchings, holders in different units,
           hub-to-port exactness) and unit/member/router ordering rows are added
  F2MK     F2M with the router graph of the shape of the known 104s: k = n - 96 hub-like routers
           without pendant doubles and 16 port-like routers with two pendant doubles and two
           hub-like neighbours (a subdivided 4-regular multigraph on k vertices)
Model(..., router_sym=False) drops the router-relabelling rows (used when the router graph is
fixed by assumptions: fixedH_control.py, plan_c_105.py).

Usage:
  python twoblock_sat.py N [--family=general|F2|F2M|F2MK] [--cap=SECONDS] [--out=runs/NAME]
  python twoblock_sat.py validate EDGEFILE          # known 104 under assumptions, must be SAT
Every SAT answer is decoded to an edge list, BFS-checked here, and written to runs/<n>_4_4_*.edges.
"""
import collections
import hashlib
import itertools
import json
import pathlib
import sys
import time

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Cadical153

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import verify  # noqa: E402


def block_vertex(b, kind, i, j=None):
    """b in {0,1}; kind 'L','R' hubs i in 0..3; 'P' port (i, j)."""
    base = 24 * b
    if kind == "L":
        return base + i
    if kind == "R":
        return base + 4 + i
    return base + 8 + 4 * i + j


class Model:
    def __init__(self, n, family="general", router_sym=True):
        self.n = n
        self.m = n - 48
        self.family = family
        self.router_sym = router_sym  # False when the router graph is fixed by assumptions
        self.pool = IDPool()
        self.solver = Cadical153()
        self.nclauses = 0
        self.W = list(range(48, n))
        self.ports = []  # (vertex id, block, l, r)
        self.hub_of_port = {}
        self.ports_of_hub = collections.defaultdict(list)
        self.block_of = {}
        self.fixed_adj = collections.defaultdict(set)
        for b in range(2):
            for i in range(4):
                self.block_of[block_vertex(b, "L", i)] = b
                self.block_of[block_vertex(b, "R", i)] = b
            for l in range(4):
                for r in range(4):
                    q = block_vertex(b, "P", l, r)
                    self.block_of[q] = b
                    hl, hr = block_vertex(b, "L", l), block_vertex(b, "R", r)
                    self.ports.append((q, b, l, r))
                    self.hub_of_port[q] = (hl, hr)
                    self.ports_of_hub[hl].append(q)
                    self.ports_of_hub[hr].append(q)
                    self.fixed_adj[q] |= {hl, hr}
                    self.fixed_adj[hl].add(q)
                    self.fixed_adj[hr].add(q)
        self.port_ids = [p[0] for p in self.ports]
        self.hub_ids = sorted(self.ports_of_hub)
        self.is_port = {q: True for q in self.port_ids}
        # fixed distances inside the blocks (BFS on block edges)
        self.fixed_dist = {}
        for s in list(self.block_of):
            dist = {s: 0}
            frontier = [s]
            while frontier:
                nxt = []
                for u in frontier:
                    for v in self.fixed_adj[u]:
                        if v not in dist:
                            dist[v] = dist[u] + 1
                            nxt.append(v)
                frontier = nxt
            for t, d in dist.items():
                self.fixed_dist[(s, t)] = d
        # star-mates: ports sharing a hub
        self.star = {q: set() for q in self.port_ids}
        for h, qs in self.ports_of_hub.items():
            for q in qs:
                self.star[q] |= set(qs) - {q}
        self.A = {}
        for c in self.W:
            for q in self.port_ids:
                self.A[(c, q)] = self.pool.id(("A", c, q))
        self.E = {}
        for u, v in itertools.combinations(self.W, 2):
            self.E[(u, v)] = self.pool.id(("E", u, v))
        self.n2cache = {}
        self.aux_count = 0
        self.build()

    # ---- helpers -------------------------------------------------------------------------------
    def add(self, clause):
        self.solver.add_clause(clause)
        self.nclauses += 1

    def evar(self, u, v):
        return self.E[(u, v)] if u < v else self.E[(v, u)]

    def adj(self, u, v):
        """True (fixed edge), False (impossible), or a variable literal."""
        if u == v:
            return False
        uw, vw = u >= 48, v >= 48
        if uw and vw:
            return self.evar(u, v)
        if uw != vw:
            c, x = (u, v) if uw else (v, u)
            if x in self.is_port:
                return self.A[(c, x)]
            return False  # W-hub
        return v in self.fixed_adj[u]

    def and_aux(self, a, b):
        """Literal equivalent (one-directionally) to a and b; a, b are literals or True."""
        if a is True:
            return b
        if b is True:
            return a
        key = ("AND", min(a, b), max(a, b))
        new = key not in self.pool.obj2id
        v = self.pool.id(key)
        if new:
            self.add([-v, a])
            self.add([-v, b])
            self.aux_count += 1
        return v

    def n2(self, u, v):
        """Literal/constant for 'walk of length <= 2 between u and v'."""
        if u == v:
            return True
        key = (min(u, v), max(u, v))
        if key in self.n2cache:
            return self.n2cache[key]
        uw, vw = u >= 48, v >= 48
        res = None
        if not uw and not vw:
            # both block vertices
            if self.block_of[u] == self.block_of[v]:
                res = self.fixed_dist[(u, v)] <= 2
            else:
                # cross block: only two ports can be at distance 2 (common holder)
                if u in self.is_port and v in self.is_port:
                    var = self.pool.id(("N2",) + key)
                    wit = [self.and_aux(self.A[(c, u)], self.A[(c, v)]) for c in self.W]
                    self.add([-var] + wit)
                    res = var
                else:
                    res = False
        else:
            if uw and vw:
                var = self.pool.id(("N2",) + key)
                wit = [self.evar(u, v)]
                for w in self.W:
                    if w != u and w != v:
                        wit.append(self.and_aux(self.evar(u, w), self.evar(w, v)))
                for q in self.port_ids:
                    wit.append(self.and_aux(self.A[(u, q)], self.A[(v, q)]))
                self.add([-var] + wit)
                res = var
            else:
                c, x = (u, v) if uw else (v, u)
                if x in self.is_port:
                    var = self.pool.id(("N2",) + key)
                    wit = [self.A[(c, x)]]
                    for w in self.W:
                        if w != c:
                            wit.append(self.and_aux(self.evar(c, w), self.A[(w, x)]))
                    self.add([-var] + wit)
                    res = var
                else:
                    # hub: c holds a port of the hub
                    var = self.pool.id(("U", c, x))
                    self.add([-var] + [self.A[(c, q)] for q in self.ports_of_hub[x]])
                    res = var
        self.n2cache[key] = res
        return res

    def within4_clause(self, u, v):
        wit = []
        for w in range(self.n):
            a = self.n2(u, w)
            if a is False:
                continue
            b = self.n2(w, v)
            if b is False:
                continue
            if a is True and b is True:
                return  # automatically within 4
            wit.append(self.and_aux(a, b))
        assert wit, (u, v)
        self.add(wit)

    def card_equal(self, lits, k):
        enc = CardEnc.equals(lits=lits, bound=k, vpool=self.pool, encoding=EncType.seqcounter)
        for cl in enc.clauses:
            self.add(cl)

    # ---- the model -----------------------------------------------------------------------------
    def build(self):
        t0 = time.time()
        # degrees
        for c in self.W:
            lits = [self.A[(c, q)] for q in self.port_ids] + [self.evar(c, v) for v in self.W if v != c]
            self.card_equal(lits, 4)
        for q in self.port_ids:
            self.card_equal([self.A[(c, q)] for c in self.W], 2)
        # family restrictions
        if self.family in ("F2", "F2M", "F2MK"):
            D = self.W[:32]
            R0 = self.W[32:]
            for c in D:
                self.card_equal([self.A[(c, q)] for q in self.port_ids], 2)
            for r in R0:
                for q in self.port_ids:
                    self.add([-self.A[(r, q)]])
            if self.family == "F2MK":
                # R0 = k hub-like routers (no pendant double) + 16 port-like routers (two pendant
                # doubles, two edges to hub-like routers): a subdivided 4-regular multigraph on k
                # vertices with 16 subdivided edges replaces the third block; k = n - 96
                k = self.n - 96
                HUBR, PORTR = R0[:k], R0[k:]
                assert len(PORTR) == 16, (k, len(PORTR))
                for h in HUBR:
                    for d in D:
                        self.add([-self.evar(d, h)])
                for p in PORTR:
                    self.card_equal([self.evar(d, p) for d in D], 2)
                    self.card_equal([self.evar(h, p) for h in HUBR], 2)
                    for p2 in PORTR:
                        if p2 > p:
                            self.add([-self.evar(p, p2)])
                # symmetry: the hub-like neighbours of PORTR[0] are HUBR[0], HUBR[1]
                if self.router_sym:
                    self.add([self.evar(HUBR[0], PORTR[0])])
                    self.add([self.evar(HUBR[1], PORTR[0])])
            if self.family in ("F2M", "F2MK"):
                for k in range(16):
                    c, c2 = D[2 * k], D[2 * k + 1]
                    self.add([self.evar(c, c2)])
                for c, c2 in itertools.combinations(D, 2):
                    if (c - 48) // 2 != (c2 - 48) // 2:
                        self.add([-self.evar(c, c2)])
                # redundant, forced by B <= 256 (RESULTS.md section 2, corollary F2M): every unit holds
                # two ports per block, they are disjoint K4,4 edges, and the two holders of a port
                # belong to different units
                for k in range(16):
                    c, c2 = D[2 * k], D[2 * k + 1]
                    for b in range(2):
                        lits = [self.A[(x, q)] for x in (c, c2) for q in self.port_ids if self.block_of[q] == b]
                        self.card_equal(lits, 2)
                    for q in self.port_ids:
                        self.add([-self.A[(c, q)], -self.A[(c2, q)]])
                        for q2 in self.star[q]:
                            if q2 > q:
                                for x in (c, c2):
                                    for y in (c, c2):
                                        self.add([-self.A[(x, q)], -self.A[(y, q2)]])
                # hub-to-port exactness (B = 256 forced): every hub sees every foreign port through
                # exactly one unit (TEMPLATE_CALCULUS 13.4 "tight mode", here at unit level)
                UA, UQ = {}, {}
                for k in range(16):
                    c, c2 = D[2 * k], D[2 * k + 1]
                    for a in self.hub_ids:
                        v = self.pool.id(("UA", k, a))
                        lits = [self.A[(x, p)] for x in (c, c2) for p in self.ports_of_hub[a]]
                        self.add([-v] + lits)
                        for l in lits:
                            self.add([-l, v])
                        UA[(k, a)] = v
                    for q in self.port_ids:
                        v = self.pool.id(("UQ", k, q))
                        lits = [self.A[(c, q)], self.A[(c2, q)]]
                        self.add([-v] + lits)
                        for l in lits:
                            self.add([-l, v])
                        UQ[(k, q)] = v
                for a in self.hub_ids:
                    for q in self.port_ids:
                        if self.block_of[q] == self.block_of[a]:
                            continue
                        prods = []
                        for k in range(16):
                            v = self.pool.id(("UAQ", k, a, q))
                            self.add([-v, UA[(k, a)]])
                            self.add([-v, UQ[(k, q)]])
                            self.add([v, -UA[(k, a)], -UQ[(k, q)]])
                            prods.append(v)
                        self.add(prods)
                        for v1, v2 in itertools.combinations(prods, 2):
                            self.add([-v1, -v2])
                # symmetry: units ordered by their smallest block-1 port (non-decreasing), and the
                # first member of a unit holds that smallest port (unit 0: W[0] holds port (0,0))
                P1 = [q for q in self.port_ids if self.block_of[q] == 0]
                MP = {}
                for k in range(16):
                    for i, q in enumerate(P1):
                        v = self.pool.id(("MP", k, q))
                        MP[(k, q)] = v
                        self.add([-v, UQ[(k, q)]])
                        for q2 in P1[:i]:
                            self.add([-v, -UQ[(k, q2)]])
                        self.add([v, -UQ[(k, q)]] + [UQ[(k, q2)] for q2 in P1[:i]])
                        self.add([-v, self.A[(D[2 * k], q)]])
                for k in range(15):
                    for i, q in enumerate(P1):
                        for q2 in P1[:i]:
                            self.add([-MP[(k, q)], -MP[(k + 1, q2)]])
        # symmetry breaking: the two holders of block-1 port (0,0) are W[0], W[1]
        q00 = block_vertex(0, "P", 0, 0)
        if self.family in ("F2", "F2M", "F2MK"):
            self.add([self.A[(self.W[0], q00)]])
            if self.family in ("F2M", "F2MK") and self.router_sym:
                # routers (port-like routers in F2MK) are interchangeable: the pendant router of
                # W[0] is the first of them, that of W[1] the first or second
                R0 = self.W[32:] if self.family == "F2M" else self.W[32 + self.n - 96:]
                self.add([self.evar(self.W[0], R0[0])])
                self.add([self.evar(self.W[1], R0[0]), self.evar(self.W[1], R0[1])])
        else:
            self.add([self.A[(self.W[0], q00)]])
            self.add([self.A[(self.W[1], q00)]])
        # diameter: every pair not inside one block
        pairs = 0
        for u in range(self.n):
            for v in range(u + 1, self.n):
                if u < 48 and v < 48 and self.block_of[u] == self.block_of[v]:
                    continue
                self.within4_clause(u, v)
                pairs += 1
        self.build_time = time.time() - t0
        self.pairs = pairs
        print(f"model n={self.n} family={self.family}: vars={self.pool.top} clauses={self.nclauses} "
              f"aux={self.aux_count} pairs={pairs} build={self.build_time:.1f}s", file=sys.stderr, flush=True)

    # ---- solving and decoding ------------------------------------------------------------------
    def solve(self, cap=None, assumptions=(), chunk=20000):
        """cap in seconds, enforced through conflict-budget chunks (the solver does not release
        the GIL, so a timer thread cannot interrupt it); None means the cap was reached."""
        t0 = time.time()
        if cap:
            res = None
            self.conflicts = 0
            while True:
                self.solver.conf_budget(chunk)
                res = self.solver.solve_limited(assumptions=list(assumptions))
                self.conflicts += chunk
                if res is not None or time.time() - t0 > cap:
                    break
        else:
            res = self.solver.solve(assumptions=list(assumptions))
        self.solve_time = time.time() - t0
        return res

    def decode(self):
        model = set(l for l in self.solver.get_model() if l > 0)
        edges = set()
        for u in self.block_of:
            for v in self.fixed_adj[u]:
                if u < v:
                    edges.add((u, v))
        for (c, q), var in self.A.items():
            if var in model:
                edges.add((min(c, q), max(c, q)))
        for (u, v), var in self.E.items():
            if var in model:
                edges.add((u, v))
        return sorted(edges)

    def assumptions_from_graph(self, edges):
        """Assumption literals fixing every A and E variable to the given graph (vertex ids in model
        indexing)."""
        es = set((min(u, v), max(u, v)) for u, v in edges)
        lits = []
        for (c, q), var in self.A.items():
            lits.append(var if (min(c, q), max(c, q)) in es else -var)
        for (u, v), var in self.E.items():
            lits.append(var if (u, v) in es else -var)
        return lits


def bfs_check(edges, n):
    adj = collections.defaultdict(set)
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    degs = collections.Counter(len(adj[v]) for v in range(n))
    diam, bad = 0, 0
    for s in range(n):
        dist = {s: 0}
        fr = [s]
        while fr:
            nx = []
            for u in fr:
                for v in adj[u]:
                    if v not in dist:
                        dist[v] = dist[u] + 1
                        nx.append(v)
            fr = nx
        if len(dist) != n:
            return {"connected": False}
        diam = max(diam, max(dist.values()))
        bad += sum(1 for t, d in dist.items() if d > 4 and t > s)
    return {"connected": True, "degrees": dict(degs), "diameter": diam, "pairs_beyond_4": bad}


def relabel_known(path, family="general"):
    """Map a known 104 (edge file) onto the model's vertex ids using two of its blocks."""
    edges = verify.load_edges_file(path)
    verts, adj, _ = verify.build(edges)
    blocks = verify.k44_blocks(verts, adj)
    keys = list(blocks)
    b1, b2 = [k for k in keys if not (k & keys[0])][0], keys[0]
    b1, b2 = keys[0], [k for k in keys if not (k & keys[0])][0]
    mapping = {}
    for b, key in enumerate((b1, b2)):
        L, R, pmap = blocks[key]
        for i, h in enumerate(L):
            mapping[h] = block_vertex(b, "L", i)
        for i, h in enumerate(R):
            mapping[h] = block_vertex(b, "R", i)
        for (l, r), q in pmap.items():
            mapping[q] = block_vertex(b, "P", L.index(l), R.index(r))
    rest = [v for v in verts if v not in mapping]
    # W[0] must hold block-1 port (0,0) (symmetry breaking), W[1] its other holder
    q00 = block_vertex(0, "P", 0, 0)
    inv = {v: k for k, v in mapping.items()}
    holders = [v for v in adj[inv[q00]] if v in rest]
    if family in ("F2M", "F2MK"):
        # doubles first, as partner pairs (2k, 2k+1), W[0] a holder of q00; routers after
        restset = set(rest)
        portset = set(inv[q] for q in [p[0] for p in [(block_vertex(b, "P", l, r), b, l, r) for b in range(2) for l in range(4) for r in range(4)]])
        doubles = [v for v in rest if adj[v] & portset]
        routers = [v for v in rest if not (adj[v] & portset)]
        partner = {}
        for c in doubles:
            nb = [u for u in adj[c] if u in restset and u in set(doubles)]
            assert len(nb) == 1, "F2M relabelling needs a perfect matching of doubles"
            partner[c] = nb[0]
        # units ordered by their smallest block-1 port (model index), first member holding it
        p1set = {inv[block_vertex(0, "P", l, r)]: block_vertex(0, "P", l, r) for l in range(4) for r in range(4)}

        def unit_key(c):
            ports = [p1set[u] for x in (c, partner[c]) for u in adj[x] if u in p1set]
            return min(ports)

        units = []
        seen = set()
        for c in doubles:
            if c in seen:
                continue
            seen |= {c, partner[c]}
            k = unit_key(c)
            first = c if any(p1set.get(u) == k for u in adj[c]) else partner[c]
            units.append((k, first, partner[first]))
        units.sort(key=lambda t: (t[0], t[1] != holders[0] and t[2] != holders[0]))
        order = []
        for k, a_, b_ in units:
            order += [a_, b_]
        rset = set(routers)
        first = [u for u in adj[order[0]] if u in rset] + [u for u in adj[order[1]] if u in rset]
        rorder = []
        for u in first + routers:
            if u not in rorder:
                rorder.append(u)
        if family == "F2MK":
            dset = set(doubles)
            portlike = [u for u in rorder if adj[u] & dset]
            hublike = [u for u in rorder if not (adj[u] & dset)]
            hfirst = [u for u in adj[portlike[0]] if u in set(hublike)]
            hublike = hfirst + [u for u in hublike if u not in hfirst]
            rorder = hublike + portlike
        order += rorder
    else:
        order = holders + [v for v in rest if v not in holders]
    for i, v in enumerate(order):
        mapping[v] = 48 + i
    return [(mapping[u], mapping[v]) for u, v in edges]


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
    if args[0] == "validate":
        edges = relabel_known(args[1], opts.get("family", "general"))
        n = len({v for e in edges for v in e})
        model = Model(n, opts.get("family", "general"))
        res = model.solve(assumptions=model.assumptions_from_graph(edges))
        print(json.dumps({"validate": args[1], "n": n, "family": model.family, "sat": res,
                          "solve_time": model.solve_time, "vars": model.pool.top, "clauses": model.nclauses}))
        return
    n = int(args[0])
    family = opts.get("family", "general")
    cap = float(opts.get("cap", 600))
    out = opts.get("out", f"runs/sat_{n}_{family}")
    model = Model(n, family)
    res = model.solve(cap=cap)
    rec = {"n": n, "family": family, "cap": cap, "result": {True: "SAT", False: "UNSAT", None: "TIMEOUT"}[res],
           "solve_time": model.solve_time, "build_time": model.build_time, "vars": model.pool.top,
           "clauses": model.nclauses, "pairs": model.pairs, "conflict_budget_used": getattr(model, "conflicts", None)}
    if res:
        edges = model.decode()
        chk = bfs_check(edges, n)
        rec["bfs"] = chk
        path = HERE / f"{out}.edges"
        path.write_text("".join(f"{u} {v}\n" for u, v in edges))
        rec["edges_file"] = str(path)
        rec["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    print(json.dumps(rec, indent=1))
    (HERE / f"{out}.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
