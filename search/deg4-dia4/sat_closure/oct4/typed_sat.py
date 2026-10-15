"""Typed exactness SAT model for block designs with a fixed connector-graph skeleton.

A skeleton lists the connector-graph components allowed by the structure theorem (potential.py /
structure_ilp.py): units (two adjacent triples whose profiles add to a permutation of
(2,2,1,1)), alternating cycles of cross doubles (adjacent doubles on complementary block pairs) and
paths t - s - t (a same-block double between two straight triples missing its block). Types inside
a component are free: the model only fixes the component shapes and the connector edges.

Variables A[c][q] (connector c holds port q). Constraints: port counts per connector (3 for
triples, 2 for doubles), every port holds exactly two connectors, the component profile
constraints, and hub-to-port exactness: for every hub a and every port q of another block exactly
one pair (c, c'') with c a slot connector of a, c'' in N[c] and c'' holding q (the routes
a-p-c-q and a-p-c-c''-q counted with multiplicity). A connector holds at most one port of a hub
(otherwise it would carry its foreign ports twice to that hub). Optional: hub-to-hub coverage
(a-p-c-q-b is the only path shape between hubs of different blocks).

Usage:
  typed_sat.py oct4 all [--nohubhub] [--atleast]     all 11 Oct x4 skeletons (units/cycles/paths)
  typed_sat.py validate104 FILE.edges                 fix a known 104 graph's own skeleton and
                                                      assignment: the exact model must be SAT
"""
import collections
import itertools
import json
import sys
import time

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Cadical153

OCT_EDGES = [(i, j) for i in range(6) for j in range(i + 1, 6) if (i, j) not in ((0, 3), (1, 4), (2, 5))]
K44_EDGES = [(i, j) for i in range(4) for j in range(4, 8)]
BASES = {"oct": (6, OCT_EDGES), "k44": (8, K44_EDGES)}


class Skeleton:
    """components: list of ("unit", 2 ids) | ("cycle", ids in cyclic order) | ("path", [t, s, t'])."""

    def __init__(self, components):
        self.components = components
        self.nc = sum(len(ids) for _, ids in components)
        assert sorted(itertools.chain.from_iterable(ids for _, ids in components)) == list(range(self.nc))
        self.edges = set()
        self.size = {}
        for kind, ids in components:
            if kind == "unit":
                (c1, c2) = ids
                self.edges.add((min(c1, c2), max(c1, c2)))
                self.size[c1] = self.size[c2] = 3
            elif kind == "cycle":
                assert len(ids) >= 4 and len(ids) % 2 == 0
                for k in range(len(ids)):
                    u, w = ids[k], ids[(k + 1) % len(ids)]
                    self.edges.add((min(u, w), max(u, w)))
                    self.size[u] = 2
            elif kind == "path":
                t, s, t2 = ids
                self.edges.add((min(t, s), max(t, s)))
                self.edges.add((min(s, t2), max(s, t2)))
                self.size[t] = self.size[t2] = 3
                self.size[s] = 2
            else:
                raise ValueError(kind)
        self.nbr = {c: set() for c in range(self.nc)}
        for u, w in self.edges:
            self.nbr[u].add(w)
            self.nbr[w].add(u)
        for c in range(self.nc):
            assert self.size[c] + len(self.nbr[c]) == 4, (c, self.size[c], self.nbr[c])


class Model:
    def __init__(self, base, nblocks, skel, exact=True, hubhub=True, aggregate=False, symbreak=False):
        nh, bedges = BASES[base]
        self.base, self.nblocks, self.skel = base, nblocks, skel
        self.pool = IDPool()
        self.cl = []
        v = self.pool.id
        self.ports = [(b, k) for b in range(nblocks) for k in range(len(bedges))]
        self.hubs = [(b, h) for b in range(nblocks) for h in range(nh)]
        pidx = {p: i for i, p in enumerate(self.ports)}
        self.pidx = pidx
        self.ports_of_hub = {hb: [pidx[(b, k)] for (b, k) in self.ports if b == hb[0] and hb[1] in bedges[k]] for hb in self.hubs}
        self.ports_of_block = {b: [pidx[p] for p in self.ports if p[0] == b] for b in range(nblocks)}
        P, NC = len(self.ports), skel.nc
        self.A = [[v(("A", c, i)) for i in range(P)] for c in range(NC)]
        for c in range(NC):
            self.eq([self.A[c][i] for i in range(P)], skel.size[c])
        for i in range(P):
            self.eq([self.A[c][i] for c in range(NC)], 2)
        # component profile constraints
        for kind, ids in skel.components:
            if kind == "unit":
                c1, c2 = ids
                for b in range(nblocks):
                    lits = [self.A[c1][i] for i in self.ports_of_block[b]] + [self.A[c2][i] for i in self.ports_of_block[b]]
                    self.cl.append(lits)  # >= 1
                    self.atmost(lits, 2)
            elif kind == "cycle":
                for k, d in enumerate(ids):
                    d2 = ids[(k + 1) % len(ids)]
                    for b in range(nblocks):
                        lits = [self.A[d][i] for i in self.ports_of_block[b]]
                        self.atmost(lits, 1)
                        # complementary neighbours: no common block
                        for i in self.ports_of_block[b]:
                            for i2 in self.ports_of_block[b]:
                                self.cl.append([-self.A[d][i], -self.A[d2][i2]])
            else:
                t, s, t2 = ids
                for b in range(nblocks):
                    for i in self.ports_of_block[b]:
                        # s: a port in block b forces a second port of s in block b
                        self.cl.append([-self.A[s][i]] + [self.A[s][i2] for i2 in self.ports_of_block[b] if i2 != i])
                        for tt in (t, t2):
                            for i2 in self.ports_of_block[b]:
                                # straight triples miss the block of s
                                self.cl.append([-self.A[s][i], -self.A[tt][i2]])
                    for tt in (t, t2):
                        self.atmost([self.A[tt][i] for i in self.ports_of_block[b]], 1)
        # U[c][a]: c holds a port of hub a (at most one)
        self.U = {}
        for c in range(NC):
            for a in self.hubs:
                lits = [self.A[c][i] for i in self.ports_of_hub[a]]
                u = v(("U", c, a))
                self.U[(c, a)] = u
                self.cl.append([-u] + lits)
                for l in lits:
                    self.cl.append([u, -l])
                for l1, l2 in itertools.combinations(lits, 2):
                    self.cl.append([-l1, -l2])
        # exactness: for hub a and foreign port q, exactly one (c, c'' in N[c]) with U[c][a] & A[c''][q]
        pairs = [(c, c2) for c in range(NC) for c2 in sorted(skel.nbr[c] | {c})]
        self.W = {}
        for a in self.hubs:
            for i, q in enumerate(self.ports):
                if q[0] == a[0]:
                    continue
                terms = []
                for c, c2 in pairs:
                    w = v(("W", c, c2, a, i))
                    self.W[(c, c2, a, i)] = w
                    self.cl += [[-w, self.U[(c, a)]], [-w, self.A[c2][i]], [w, -self.U[(c, a)], -self.A[c2][i]]]
                    terms.append(w)
                if exact:
                    self.eq(terms, 1)
                else:
                    self.cl.append(terms)
        if aggregate == "wsum" and exact:
            # the sum of the six exactly-one constraints of the hubs of block b for a fixed foreign
            # port q: exactly six route literals (a consequence of the above; added so that the
            # solver sees the double count  f_b(h1) + f_b(h2) = 3  of the holders h1, h2 of q)
            for i, q in enumerate(self.ports):
                for b in range(nblocks):
                    if b == q[0]:
                        continue
                    lits = [self.W[(c, c2, a, i)] for a in self.hubs if a[0] == b for (c, c2) in pairs]
                    self.eq(lits, 6)
        if aggregate == "fcount" and exact:
            # the same double count with explicit carried counts: X[c][b][v] <- (c has exactly v
            # ports in block b), F[h][b][w] <- (N[h] carries exactly w ports of block b), both encoded
            # in the forcing direction only (v, w <= 3; larger counts force nothing, which only
            # weakens the redundant constraint), and for every port q of another block and every pair
            # of potential holders: both hold q -> F values add to 3.
            X = {}
            for c in range(NC):
                for b in range(nblocks):
                    lits = [self.A[c][i] for i in self.ports_of_block[b]]
                    prev = {0: None}  # state after 0 literals: count 0 is certain (None = true)
                    for k, l in enumerate(lits):
                        cur = {}
                        for val in range(0, min(k + 1, 3) + 1):
                            cur[val] = v(("X", c, b, k, val))
                        for val, sv in cur.items():
                            # count stays val (literal false) or increases to val (literal true)
                            if val in prev:
                                self.cl.append(([-prev[val]] if prev[val] is not None else []) + [l, sv])
                            if val - 1 in prev:
                                self.cl.append(([-prev[val - 1]] if prev[val - 1] is not None else []) + [-l, sv])
                        prev = cur
                    for val in range(4):
                        X[(c, b, val)] = prev.get(val)
            F = {}
            for h in range(NC):
                nb = sorted(skel.nbr[h] | {h})
                for b in range(nblocks):
                    for w in range(4):
                        F[(h, b, w)] = v(("F", h, b, w))
                    for vals in itertools.product(range(4), repeat=len(nb)):
                        w = sum(vals)
                        if w > 3:
                            continue
                        if any(X[(c, b, val)] is None and val != 0 for c, val in zip(nb, vals)):
                            continue
                        self.cl.append([-X[(c, b, val)] for c, val in zip(nb, vals) if X[(c, b, val)] is not None] + [F[(h, b, w)]])
            # the holders h1, h2 of q carry 2 (f_b(h1) + f_b(h2)) routes from the nh hubs of block b
            target = nh // 2
            for i, q in enumerate(self.ports):
                for b in range(nblocks):
                    if b == q[0]:
                        continue
                    for h1 in range(NC):
                        for h2 in range(h1 + 1, NC):
                            for w1 in range(4):
                                for w2 in range(4):
                                    if w1 + w2 != target:
                                        self.cl.append([-self.A[h1][i], -self.A[h2][i], -F[(h1, b, w1)], -F[(h2, b, w2)]])
        if symbreak:
            # block permutations: unit 0 (if any) is heavy in blocks 0 and 1
            units = [ids for kind, ids in skel.components if kind == "unit"]
            if units:
                c0, c1 = units[0]
                for b in (0, 1):
                    lits = [self.A[c0][i] for i in self.ports_of_block[b]] + [self.A[c1][i] for i in self.ports_of_block[b]]
                    for cl in CardEnc.atleast(lits=lits, bound=2, vpool=self.pool, encoding=EncType.seqcounter).clauses:
                        self.cl.append(cl)
        if symbreak:
            # block automorphisms act independently and transitively on each block's ports, and the
            # connectors of a unit / the units / the paths are interchangeable in the skeleton
            units = [ids for kind, ids in skel.components if kind == "unit"]
            paths = [ids for kind, ids in skel.components if kind == "path"]
            if units:
                c0, c1 = units[0]
                # one connector of unit 0 holds a block-0 port: make it c0 and the port (0, 0)
                self.cl.append([self.A[c0][pidx[(0, 0)]]])
                # unit 0 holds a port in every block: make it port 0 of that block
                for b in range(1, nblocks):
                    self.cl.append([self.A[c0][pidx[(b, 0)]], self.A[c1][pidx[(b, 0)]]])
            elif paths:
                t, s, t2 = paths[0]
                # block permutation: the same-block double of path 0 sits in block 0, at port (0, 0);
                # its straight triple t holds one port in each other block: port 0 of that block
                self.cl.append([self.A[s][pidx[(0, 0)]]])
                for b in range(1, nblocks):
                    self.cl.append([self.A[t][pidx[(b, 0)]]])
        if hubhub:
            for a in self.hubs:
                for b in self.hubs:
                    if a >= b or a[0] == b[0]:
                        continue
                    terms = []
                    for c in range(NC):
                        w = v(("HH", c, a, b))
                        self.cl += [[-w, self.U[(c, a)]], [-w, self.U[(c, b)]]]
                        terms.append(w)
                    self.cl.append(terms)

    def eq(self, lits, k):
        for cl in CardEnc.equals(lits=lits, bound=k, vpool=self.pool, encoding=EncType.seqcounter).clauses:
            self.cl.append(cl)

    def atmost(self, lits, k):
        for cl in CardEnc.atmost(lits=lits, bound=k, vpool=self.pool, encoding=EncType.seqcounter).clauses:
            self.cl.append(cl)

    def decode(self, model):
        pos = set(l for l in model if l > 0)
        nh, bedges = BASES[self.base]
        edges = set()
        for (b, k) in self.ports:
            u, w = bedges[k]
            edges.add((("P", b, k), ("H", b, u)))
            edges.add((("P", b, k), ("H", b, w)))
        for c in range(self.skel.nc):
            for i, p in enumerate(self.ports):
                if self.A[c][i] in pos:
                    edges.add((("C", c), ("P",) + p))
        for u, w in self.skel.edges:
            edges.add((("C", u), ("C", w)))
        return sorted((str(x), str(y)) for x, y in edges)


def bfs_check(edges):
    adj = collections.defaultdict(set)
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    verts = sorted(adj)
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
        bad += len(verts) - len(dist)
        diam = max(diam, max(dist.values()))
        bad += sum(1 for u, d in dist.items() if d > 4 and u > v)
    return {"n": len(verts), "degrees": dict(collections.Counter(len(adj[v]) for v in verts)), "diameter": diam, "bad_pairs": bad}


def oct4_skeletons():
    """All connector-graph skeletons of a tight Oct x4 design: u units, 12 - u paths t-s-t (24 triples
    in total), u cycle doubles (12 doubles in total) in cycles of even length >= 4."""
    def even_partitions(n, mn=4):
        if n == 0:
            yield []
            return
        for first in range(mn, n + 1, 2):
            for rest in even_partitions(n - first, first):
                yield [first] + rest

    out = []
    for u in range(0, 13):
        k = 12 - u
        for cyc in even_partitions(u):
            comps = []
            nxt = 0
            for _ in range(u):
                comps.append(("unit", [nxt, nxt + 1]))
                nxt += 2
            for _ in range(k):
                comps.append(("path", [nxt, nxt + 1, nxt + 2]))
                nxt += 3
            for L in cyc:
                comps.append(("cycle", list(range(nxt, nxt + L))))
                nxt += L
            assert nxt == 36
            out.append((f"u{u}_k{k}_cyc{'-'.join(map(str, cyc)) or 'none'}", comps))
    return out


def run_oct4(args):
    exact = "--atleast" not in args
    hubhub = "--nohubhub" not in args
    aggregate = next((a.split("=")[1] for a in args if a.startswith("--aggregate=")), False)
    symbreak = "--symbreak" in args
    limit = next((int(a.split("=")[1]) for a in args if a.startswith("--limit=")), None)
    only = [a for a in args if not a.startswith("--") and a != "all"]
    results = []
    for name, comps in oct4_skeletons():
        if only and name not in only:
            continue
        t0 = time.time()
        m = Model("oct", 4, Skeleton(comps), exact=exact, hubhub=hubhub, aggregate=aggregate, symbreak=symbreak)
        rec = {"skeleton": name, "exact": exact, "hubhub": hubhub, "aggregate": aggregate, "symbreak": symbreak,
               "vars": m.pool.top, "clauses": len(m.cl), "build_s": round(time.time() - t0, 1)}
        s = Cadical153(bootstrap_with=m.cl)
        t1 = time.time()
        if limit:
            # conflict budget as a proxy for a time limit (python-sat exposes no wall-clock limit)
            s.conf_budget(limit)
            r = s.solve_limited()
        else:
            r = s.solve()
        rec["result"] = r
        rec["solve_s"] = round(time.time() - t1, 1)
        if r:
            edges = m.decode(s.get_model())
            rec["check"] = bfs_check(edges)
            json.dump({"edges": edges, "check": rec["check"]}, open(f"runs/typed_{name}_solution.json", "w"))
        s.delete()
        print(json.dumps(rec), file=sys.stderr, flush=True)
        results.append(rec)
    json.dump(results, sys.stdout, indent=1)


def validate104(path):
    """Decompose a known 104 graph (three K4,4 blocks), build the exact model on its own skeleton and
    check that its own port assignment satisfies it (SAT under assumptions), and that it is UNSAT
    after moving one port."""
    sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent.parent / "template_calculus"))
    import verify  # template_calculus/verify.py

    edges = verify.load_edges_file(path)
    verts, adj, problems = verify.build(edges)
    assert not problems
    blocks = verify.k44_blocks(verts, adj)
    decs = verify.decompositions(blocks, len(verts))
    assert len(decs) == 1, len(decs)
    trip = decs[0]
    # map graph vertices to model hubs/ports: block b, hubs L (0..3) then R (4..7), port (l, r) -> K44 edge index
    hub_of, port_of = {}, {}
    for b, key in enumerate(trip):
        L, R, pmap = blocks[key]
        for i, l in enumerate(L):
            hub_of[l] = (b, i)
        for i, r in enumerate(R):
            hub_of[r] = (b, 4 + i)
        for (l, r), q in pmap.items():
            port_of[q] = (b, K44_EDGES.index((L.index(l), 4 + R.index(r))))
    connectors = [x for x in verts if x not in hub_of and x not in port_of]
    cid = {c: i for i, c in enumerate(connectors)}
    cadj = {c: [u for u in adj[c] if u in cid] for c in connectors}
    assert all(len(cadj[c]) == 1 for c in connectors), "known 104 graphs are perfect matchings of triples"
    comps, seen = [], set()
    for c in connectors:
        if c in seen:
            continue
        c2 = cadj[c][0]
        seen |= {c, c2}
        comps.append(("unit", [cid[c], cid[c2]]))
    # relabel so that ids are 0..31 in component order
    order = list(itertools.chain.from_iterable(ids for _, ids in comps))
    new = {old: i for i, old in enumerate(order)}
    comps = [(k, [new[i] for i in ids]) for k, ids in comps]
    cid = {c: new[i] for c, i in cid.items()}
    skel = Skeleton(comps)
    aggregate = next((a.split("=")[1] for a in sys.argv if a.startswith("--aggregate=")), False)
    m = Model("k44", 3, skel, exact=True, hubhub=True, aggregate=aggregate, symbreak=False)
    assumptions = []
    truth = set()
    for c in connectors:
        for u in adj[c]:
            if u in port_of:
                truth.add((cid[c], m.pidx[port_of[u]]))
    if "--symbreak" in sys.argv:
        # the symmetry-breaking clauses are tested separately: without assumptions the model must stay
        # SAT (the known graph can be relabelled to satisfy them)
        m2 = Model("k44", 3, skel, exact=True, hubhub=True, aggregate=aggregate, symbreak=True)
        s2 = Cadical153(bootstrap_with=m2.cl)
        symbreak_sat = s2.solve()
        s2.delete()
    else:
        symbreak_sat = None
    for c in range(skel.nc):
        for i in range(len(m.ports)):
            assumptions.append(m.A[c][i] if (c, i) in truth else -m.A[c][i])
    s = Cadical153(bootstrap_with=m.cl)
    ok = s.solve(assumptions=assumptions)
    rec = {"file": path, "connectors": len(connectors), "aggregate": aggregate, "vars": m.pool.top, "clauses": len(m.cl),
           "own_assignment_sat": ok, "free_model_with_symbreak_sat": symbreak_sat}
    # perturbation: swap the ports of two connectors in one block (keeps counts, breaks exactness)
    (c0, i0) = next(iter(sorted(truth)))
    other = next((c, i) for (c, i) in sorted(truth) if c != c0 and m.ports[i][0] == m.ports[i0][0] and (c0, i) not in truth and (c, i0) not in truth)
    perturbed = (truth - {(c0, i0), other}) | {(c0, other[1]), (other[0], i0)}
    assumptions2 = [m.A[c][i] if (c, i) in perturbed else -m.A[c][i] for c in range(skel.nc) for i in range(len(m.ports))]
    rec["perturbed_assignment_sat"] = s.solve(assumptions=assumptions2)
    t0 = time.time()
    rec["free_sat"] = s.solve()
    rec["free_solve_s"] = round(time.time() - t0, 1)
    if rec["free_sat"]:
        rec["free_check"] = bfs_check(m.decode(s.get_model()))
    s.delete()
    print(json.dumps(rec), file=sys.stderr, flush=True)
    json.dump(rec, sys.stdout, indent=1)


if __name__ == "__main__":
    if sys.argv[1] == "oct4":
        run_oct4(sys.argv[2:])
    elif sys.argv[1] == "validate104":
        validate104(sys.argv[2])
