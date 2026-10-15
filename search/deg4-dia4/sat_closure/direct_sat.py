"""Direct SAT hunt for a diameter-4 graph with max degree 4 built from free blocks plus connectors.

Blocks: subdivisions of K5, the octahedron K_{2,2,2}, or K4,4 (hubs = base vertices, ports = base
edges; every hub is adjacent to or on every base edge, so all pairs inside a block are within 4).
Connectors: m = n - (block vertices) extra vertices; every port has two free slots that connectors
fill; connectors may also be adjacent to each other. Variables: A[c][p] (connector c holds port p),
E[c][c'] (connector edge). Constraints: every port has exactly two connectors, every connector has
degree exactly four, and every pair of vertices not inside one block is within distance 4, encoded
through the path shapes available in this architecture (hubs touch only ports). Any model is decoded
and verified by BFS before being reported. This is a satisfiability hunt, not an impossibility proof.
"""
import collections
import itertools
import json
import sys
import time

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Cadical153

BASES = {
    "K5": [(i, j) for i in range(5) for j in range(i + 1, 5)],
    "Oct": [(i, j) for i in range(6) for j in range(i + 1, 6) if j != i + 3 - (3 if i >= 3 else 0) and (i, j) not in ((0, 3), (1, 4), (2, 5))],
    "K44": [(i, 4 + j) for i in range(4) for j in range(4)],
}
BASE_N = {"K5": 5, "Oct": 6, "K44": 8}


def check_free(kind):
    """Verify the subdivided base has all internal pairs within 4 and hubs of degree 4."""
    edges = BASES[kind]
    n = BASE_N[kind]
    adj = collections.defaultdict(set)
    for k, (u, v) in enumerate(edges):
        port = ("p", k)
        adj[("h", u)].add(port)
        adj[("h", v)].add(port)
        adj[port] |= {("h", u), ("h", v)}
    assert all(len(adj[("h", u)]) == 4 for u in range(n)), kind
    verts = list(adj)
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
        assert len(dist) == len(verts) and max(dist.values()) <= 4, (kind, v, max(dist.values()))
    return len(edges)


class Direct:
    def __init__(self, blocks, n):
        self.blocks = blocks
        self.pool = IDPool()
        self.cnf = []
        self.hubs, self.ports = [], []
        self.hubs_of_port = {}
        self.ports_of_hub = collections.defaultdict(list)
        self.block_of = {}
        self.base_nbrs = collections.defaultdict(set)
        for bi, kind in enumerate(blocks):
            check_free(kind)
            for u in range(BASE_N[kind]):
                h = ("H", bi, u)
                self.hubs.append(h)
                self.block_of[h] = bi
            for k, (u, v) in enumerate(BASES[kind]):
                p = ("P", bi, k)
                self.ports.append(p)
                self.block_of[p] = bi
                self.hubs_of_port[p] = (("H", bi, u), ("H", bi, v))
                self.ports_of_hub[("H", bi, u)].append(p)
                self.ports_of_hub[("H", bi, v)].append(p)
                self.base_nbrs[("H", bi, u)].add(("H", bi, v))
                self.base_nbrs[("H", bi, v)].add(("H", bi, u))
        self.m = n - len(self.hubs) - len(self.ports)
        assert self.m > 0
        self.P = len(self.ports)
        self.conns = list(range(self.m))
        v = self.pool.id
        self.A = [[v(("A", c, i)) for i in range(self.P)] for c in self.conns]
        self.E = {}
        for c in self.conns:
            for c2 in self.conns:
                if c < c2:
                    self.E[(c, c2)] = self.E[(c2, c)] = v(("E", c, c2))
        self.T = v(("TRUE",))
        self.cnf.append([self.T])
        self.build()

    def eq(self, lits, k):
        for cl in CardEnc.equals(lits=lits, bound=k, vpool=self.pool, encoding=EncType.seqcounter).clauses:
            self.cnf.append(cl)

    def AND(self, *lits):
        lits = tuple(sorted(set(l for l in lits if l != self.T)))
        if not lits:
            return self.T
        if len(lits) == 1:
            return lits[0]
        key = ("and",) + lits
        w = self.pool.id(key)
        if key not in self._defined:
            self._defined.add(key)
            for l in lits:
                self.cnf.append([-w, l])
            self.cnf.append([w] + [-l for l in lits])
        return w

    def OR(self, lits):
        lits = [l for l in set(lits) if l is not None]
        if self.T in lits:
            return self.T
        if not lits:
            return None
        if len(lits) == 1:
            return lits[0]
        key = ("or",) + tuple(sorted(lits))
        w = self.pool.id(key)
        if key not in self._defined:
            self._defined.add(key)
            self.cnf.append([-w] + lits)
            for l in lits:
                self.cnf.append([w, -l])
        return w

    def build(self):
        self._defined = set()
        A, E, C = self.A, self.E, self.conns
        pidx = {p: i for i, p in enumerate(self.ports)}
        hidx = {h: i for i, h in enumerate(self.hubs)}
        # degrees
        for i in range(self.P):
            self.eq([A[c][i] for c in C], 2)
        for c in C:
            self.eq([A[c][i] for i in range(self.P)] + [E[(c, c2)] for c2 in C if c2 != c], 4)
        # derived: U[c][h] = c holds a port of hub h
        U = [[self.OR([A[c][pidx[p]] for p in self.ports_of_hub[h]]) for h in self.hubs] for c in C]
        # S[c][c2] = c and c2 share a port
        S = {}
        for c in C:
            for c2 in C:
                if c < c2:
                    S[(c, c2)] = S[(c2, c)] = self.OR([self.AND(A[c][i], A[c2][i]) for i in range(self.P)])
        # E2[c][c2] = connector path of length 2 (c != c2)
        E2 = {}
        for c in C:
            for c2 in C:
                if c < c2:
                    E2[(c, c2)] = E2[(c2, c)] = self.OR([self.AND(E[(c, x)], E[(x, c2)]) for x in C if x not in (c, c2)])
        # Tcarry[c][q] = c holds q or a connector neighbour of c holds q
        Tc = [[self.OR([A[c][i]] + [self.AND(E[(c, c2)], A[c2][i]) for c2 in C if c2 != c]) for i in range(self.P)] for c in C]
        # T2[c][q] = some connector within connector-distance 2 of c holds q (via E2) or Tc
        T2 = [[self.OR([Tc[c][i]] + [self.AND(E2[(c, c2)], A[c2][i]) for c2 in C if c2 != c]) for i in range(self.P)] for c in C]
        # Vs[c][q] = a port-mate of c holds q
        Vs = [[self.OR([self.AND(S[(c, c2)], A[c2][i]) for c2 in C if c2 != c]) for i in range(self.P)] for c in C]
        n_pairs = 0
        # hub-hub across blocks: a-p-c-q-b
        for a in self.hubs:
            for b in self.hubs:
                if hidx[a] >= hidx[b] or self.block_of[a] == self.block_of[b]:
                    continue
                lits = [self.AND(A[c][pidx[p]], A[c][pidx[q]]) for c in C for p in self.ports_of_hub[a] for q in self.ports_of_hub[b]]
                self.cnf.append([l for l in lits])
                n_pairs += 1
        # hub-port across blocks: a-p-c-q, a-p-c-c'-q
        for a in self.hubs:
            for q in self.ports:
                if self.block_of[a] == self.block_of[q]:
                    continue
                lits = [self.AND(A[c][pidx[p]], Tc[c][pidx[q]]) for c in C for p in self.ports_of_hub[a]]
                self.cnf.append(lits)
                n_pairs += 1
        # port-port across blocks
        for p in self.ports:
            for q in self.ports:
                if pidx[p] >= pidx[q] or self.block_of[p] == self.block_of[q]:
                    continue
                ip, iq = pidx[p], pidx[q]
                lits = []
                for c in C:
                    lits.append(self.AND(A[c][ip], T2[c][iq]))          # p-c-q, p-c-c'-q, p-c-c'-c''-q
                    lits.append(self.AND(A[c][ip], Vs[c][iq]))          # p-c-y-c'-q
                    for b in self.hubs_of_port[q]:
                        lits.append(self.AND(A[c][ip], U[c][hidx[b]]))  # p-c-y-b-q
                    for a in self.hubs_of_port[p]:
                        lits.append(self.AND(U[c][hidx[a]], A[c][iq]))  # p-a-p'-c-q
                self.cnf.append(lits)
                n_pairs += 1
        # hub-connector: a-p-c; a-p-c'-c; a-p-c'-c''-c; a-p-c'-q-c; a-p-b-p'-c (b base-neighbour of a)
        for a in self.hubs:
            ia = hidx[a]
            nb = [hidx[b] for b in self.base_nbrs[a]]
            for c in C:
                lits = [U[c][ia]] + [U[c][ib] for ib in nb]
                for c2 in C:
                    if c2 == c:
                        continue
                    lits.append(self.AND(U[c2][ia], E[(c2, c)]))
                    lits.append(self.AND(U[c2][ia], E2[(c2, c)]))
                    lits.append(self.AND(U[c2][ia], S[(c2, c)]))
                self.cnf.append(lits)
                n_pairs += 1
        # port-connector
        for p in self.ports:
            ip = pidx[p]
            hs = [hidx[h] for h in self.hubs_of_port[p]]
            for c in C:
                lits = [A[c][ip]] + [U[c][ih] for ih in hs]
                for c2 in C:
                    if c2 == c:
                        continue
                    lits.append(self.AND(A[c2][ip], E[(c2, c)]))           # p-c'-c
                    lits.append(self.AND(A[c2][ip], E2[(c2, c)]))          # p-c'-c''-c
                    lits.append(self.AND(A[c2][ip], S[(c2, c)]))           # p-c'-q-c
                    for ih in hs:
                        lits.append(self.AND(U[c2][ih], E[(c2, c)]))       # p-a-p'-c'-c
                    for c3 in C:
                        if c3 in (c, c2):
                            continue
                        lits.append(self.AND(A[c2][ip], E[(c2, c3)], E2[(c3, c)]))   # length 4 connector path
                        lits.append(self.AND(A[c2][ip], S[(c2, c3)], E[(c3, c)]))    # p-c'-q-c''-c
                        lits.append(self.AND(A[c2][ip], E[(c2, c3)], S[(c3, c)]))    # p-c'-c''-q-c
                self.cnf.append(lits)
                n_pairs += 1
        # connector-connector
        for c in C:
            for c2 in C:
                if c >= c2:
                    continue
                lits = [E[(c, c2)], E2[(c, c2)], S[(c, c2)]]
                lits += [self.AND(U[c][ih], U[c2][ih]) for ih in range(len(self.hubs))]  # c-p-a-p'-c2
                for c3 in C:
                    if c3 in (c, c2):
                        continue
                    lits.append(self.AND(E[(c, c3)], E2[(c3, c2)]))      # length 3
                    lits.append(self.AND(S[(c, c3)], E[(c3, c2)]))       # c-q-c3-c2
                    lits.append(self.AND(E[(c, c3)], S[(c3, c2)]))
                    lits.append(self.AND(S[(c, c3)], S[(c3, c2)]))       # c-q-c3-q'-c2
                    lits.append(self.AND(E2[(c, c3)], E2[(c3, c2)]))     # length 4 connector path
                    lits.append(self.AND(E[(c, c3)], S[(c3, c2)]))
                    lits.append(self.AND(S[(c, c3)], E2[(c3, c2)]))      # c-q-c3-x-c2
                    lits.append(self.AND(E2[(c, c3)], S[(c3, c2)]))
                self.cnf.append(lits)
                n_pairs += 1
        self.n_pairs = n_pairs

    def decode(self, model):
        pos = set(l for l in model if l > 0)
        edges = set()
        for p in self.ports:
            for h in self.hubs_of_port[p]:
                edges.add((str(p), str(h)))
        for c in self.conns:
            for i, p in enumerate(self.ports):
                if self.A[c][i] in pos:
                    edges.add((str(("C", c)), str(p)))
            for c2 in self.conns:
                if c < c2 and self.E[(c, c2)] in pos:
                    edges.add((str(("C", c)), str(("C", c2))))
        return sorted(edges)


def bfs_check(edges):
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
        bad += len(verts) - len(dist)
        diam = max(diam, max(dist.values()))
        bad += sum(1 for u, d in dist.items() if d > 4 and u > v)
    return {"n": len(verts), "degrees": dict(degs), "diameter": diam, "bad_pairs": bad}


def main():
    blocks = sys.argv[1].split("+")
    n = int(sys.argv[2])
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else 600
    t0 = time.time()
    D = Direct(blocks, n)
    rec = {"blocks": blocks, "n": n, "connectors": D.m, "vars": D.pool.top, "clauses": len(D.cnf),
           "pairs": D.n_pairs, "build_s": round(time.time() - t0, 1)}
    print(json.dumps(rec), file=sys.stderr, flush=True)
    solver = Cadical153(bootstrap_with=D.cnf)
    import threading
    timer = threading.Timer(limit, lambda: solver.interrupt())
    timer.start()
    t1 = time.time()
    res = solver.solve_limited(expect_interrupt=True)
    timer.cancel()
    rec["solve_s"] = round(time.time() - t1, 1)
    rec["result"] = {True: "SAT", False: "UNSAT", None: "TIMEOUT"}[res]
    if res:
        edges = D.decode(solver.get_model())
        chk = bfs_check(edges)
        rec["verified"] = chk
        tag = "+".join(blocks)
        with open(f"runs/direct_{tag}_n{n}.json", "w") as f:
            json.dump({"blocks": blocks, "n": n, "edges": edges, "check": chk}, f)
        if chk["diameter"] == 4 and chk["bad_pairs"] == 0:
            print("!!! VERIFIED DIAMETER-4 GRAPH, n =", chk["n"], file=sys.stderr, flush=True)
    print(json.dumps(rec), file=sys.stderr, flush=True)
    json.dump(rec, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
