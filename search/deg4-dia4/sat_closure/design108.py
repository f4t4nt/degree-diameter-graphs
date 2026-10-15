"""Exactness model for 108 = four octahedron blocks + 36 connectors (profile-free).

Four subdivided octahedra (6 hubs, 12 ports, 24 slots each; 72 vertices, 96 slots) plus 36
connectors give n = 108. The hub-to-port budget relaxation has zero slack there (blocks_relax.py),
so any such diameter-4 graph is tight: every hub reaches every port of the other three blocks
exactly once along a-p-c-q or a-p-c-c'-q, and in particular no connector holds two ports of one hub
and no connector has two routes to the same port (itself and a neighbour, or two neighbours).

Variables: A[c][q] (connector c holds port q), E[c][c'] (connector edge). Constraints: every port
holds exactly two connectors; every connector has degree exactly four (ports plus connector edges);
U[c][a] = c holds a port of hub a, at most one; T[c][q] = exactly one of {A[c][q]} u {E[c][c'] & A[c'][q]};
for every hub a and every port q of another block, exactly one c with U[c][a] & T[c][q]; every pair
of hubs in different blocks shares a connector (hub-to-hub coverage). Profiles are free (routers
allowed). UNSAT closes the whole 4-octahedron family at 108; SAT solutions are decoded and checked
by BFS and would then need the remaining pair classes.
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
NB, NC = 4, 36


class Model108:
    def __init__(self, exact=True, shape=False):
        self.pool = IDPool()
        self.cl = []
        v = self.pool.id
        self.ports = [(b, k) for b in range(NB) for k in range(12)]
        self.hubs = [(b, h) for b in range(NB) for h in range(6)]
        self.ports_of_hub = {hb: [(b, k) for (b, k) in self.ports if b == hb[0] and hb[1] in OCT_EDGES[k]] for hb in self.hubs}
        P = len(self.ports)
        self.A = [[v(("A", c, i)) for i in range(P)] for c in range(NC)]
        self.E = {}
        for c in range(NC):
            for c2 in range(c + 1, NC):
                self.E[(c, c2)] = self.E[(c2, c)] = v(("E", c, c2))
        for i in range(P):
            self.eq([self.A[c][i] for c in range(NC)], 2)
        for c in range(NC):
            self.eq([self.A[c][i] for i in range(P)] + [self.E[(c, c2)] for c2 in range(NC) if c2 != c], 4)
        pidx = {p: i for i, p in enumerate(self.ports)}
        # U[c][a]: exactly-one-or-zero over the 4 ports of hub a
        self.U = {}
        for c in range(NC):
            for a in self.hubs:
                lits = [self.A[c][pidx[p]] for p in self.ports_of_hub[a]]
                u = v(("U", c, a))
                self.U[(c, a)] = u
                self.cl.append([-u] + lits)
                for l in lits:
                    self.cl.append([u, -l])
                if exact:
                    for l1, l2 in itertools.combinations(lits, 2):
                        self.cl.append([-l1, -l2])
        # T[c][q]: exactly one route from c to q (own port or via one neighbour)
        self.T = {}
        for c in range(NC):
            for i in range(P):
                terms = [self.A[c][i]]
                for c2 in range(NC):
                    if c2 == c:
                        continue
                    w = v(("W", c, c2, i))
                    self.cl += [[-w, self.E[(c, c2)]], [-w, self.A[c2][i]], [w, -self.E[(c, c2)], -self.A[c2][i]]]
                    terms.append(w)
                t = v(("T", c, i))
                self.T[(c, i)] = t
                self.cl.append([-t] + terms)
                for x in terms:
                    self.cl.append([t, -x])
                if exact:
                    for cl in CardEnc.atmost(lits=terms, bound=1, vpool=self.pool, encoding=EncType.seqcounter).clauses:
                        self.cl.append(cl)
        # hub -> foreign port exactly once
        for a in self.hubs:
            for i, q in enumerate(self.ports):
                if q[0] == a[0]:
                    continue
                terms = []
                for c in range(NC):
                    w = v(("R", c, a, i))
                    self.cl += [[-w, self.U[(c, a)]], [-w, self.T[(c, i)]], [w, -self.U[(c, a)], -self.T[(c, i)]]]
                    terms.append(w)
                if exact:
                    self.eq(terms, 1)
                else:
                    self.cl.append(terms)
        # hub-hub coverage across blocks
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
        # implied by exactness: a connector adjacent to c cannot hold a port of a hub that c serves
        # (the hub would reach c' twice), and two connectors sharing a port are not adjacent
        for c in range(NC):
            for c2 in range(NC):
                if c2 == c:
                    continue
                for a in self.hubs:
                    self.cl.append([-self.U[(c, a)], -self.E[(c, c2)], -self.U[(c2, a)]])
                if c < c2:
                    for i in range(P):
                        self.cl.append([-self.A[c][i], -self.A[c2][i], -self.E[(c, c2)]])
        # symmetry breaking (connectors are interchangeable): first held port index non-decreasing
        # along connectors, routers last. F[c][i] = c holds some port with index <= i.
        F = [[v(("F", c, i)) for i in range(P)] for c in range(NC)]
        for c in range(NC):
            for i in range(P):
                lits = [self.A[c][j] for j in range(i + 1)]
                self.cl.append([-F[c][i]] + lits)
                for l in lits:
                    self.cl.append([F[c][i], -l])
        for c in range(NC - 1):
            for i in range(P):
                self.cl.append([-F[c + 1][i], F[c][i]])
        # port 0 is held by exactly two connectors, which are then connectors 0 and 1
        self.cl.append([self.A[0][0]])
        self.cl.append([self.A[1][0]])
        if shape:
            # every tight shape has exactly 24 connectors with 3 ports and 12 with 2 ports
            # (ILP enumeration of all 49 tight type-count vectors, runs/design108_shapes.txt):
            # forbid 0, 1 and 4 ports; with 96 slots and 36 connectors the counts are then forced.
            for c in range(NC):
                for cl in CardEnc.atleast(lits=[self.A[c][i] for i in range(P)], bound=2, vpool=self.pool,
                                          encoding=EncType.seqcounter).clauses:
                    self.cl.append(cl)
                for cl in CardEnc.atmost(lits=[self.A[c][i] for i in range(P)], bound=3, vpool=self.pool,
                                         encoding=EncType.seqcounter).clauses:
                    self.cl.append(cl)

    def eq(self, lits, k):
        for cl in CardEnc.equals(lits=lits, bound=k, vpool=self.pool, encoding=EncType.seqcounter).clauses:
            self.cl.append(cl)

    def decode(self, model):
        pos = set(l for l in model if l > 0)
        edges = set()
        for (b, k) in self.ports:
            u, w = OCT_EDGES[k]
            edges.add((("P", b, k), ("H", b, u)))
            edges.add((("P", b, k), ("H", b, w)))
        for c in range(NC):
            for i, p in enumerate(self.ports):
                if self.A[c][i] in pos:
                    edges.add((("C", c), ("P",) + p))
            for c2 in range(c + 1, NC):
                if self.E[(c, c2)] in pos:
                    edges.add((("C", c), ("C", c2)))
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


def main():
    exact = (sys.argv[1] != "relaxed") if len(sys.argv) > 1 else True
    shape = len(sys.argv) > 2 and sys.argv[2] == "shape"
    t0 = time.time()
    m = Model108(exact=exact, shape=shape)
    rec = {"exact": exact, "shape": shape, "vars": m.pool.top, "clauses": len(m.cl), "build_s": round(time.time() - t0, 1)}
    print(json.dumps(rec), file=sys.stderr, flush=True)
    s = Cadical153(bootstrap_with=m.cl)
    t1 = time.time()
    r = s.solve()
    rec["result"] = r
    rec["solve_s"] = round(time.time() - t1, 1)
    if r:
        edges = m.decode(s.get_model())
        rec["check"] = bfs_check(edges)
        json.dump({"edges": edges, "check": rec["check"]}, open("runs/design108_solution.json", "w"))
    s.delete()
    print(json.dumps(rec), file=sys.stderr, flush=True)
    json.dump(rec, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
