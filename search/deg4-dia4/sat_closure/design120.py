"""Exact model for the 120-vertex all-cross-double shape of the three-K4,4-block family.

Shape (counting, zero slack): 48 doubles, 16 of each type {1,2}, {1,3}, {2,3}; every double has one
port in each of its two blocks and two connector edges to doubles of the other types (same-type
edges lose budget and are excluded). Hub-to-hub exactness forces every port to hold exactly one
double of each of its two types, so the {i,j}-doubles define a bijection phi_ij: P_i -> P_j that maps
every hub star to a perfect matching ("OA bijection"; equivalently two orthogonal Latin squares of
order 4). Zero slack forces hub-to-port exactness: every hub reaches every port of each other block
exactly once along a-p-d-q or a-p-d-d'-q.

Model: phi_12 and phi_13 fixed up to symmetry (enumerated OA bijections reduced by K4,4
automorphisms), phi_23 as permutation variables with the OA constraint, edge variables between
doubles of different types with degree exactly two, and the twelve families of exactly-one reach
constraints. Every solution is turned into the 120-vertex graph and checked by BFS for diameter 4.
"""
import collections
import itertools
import json
import sys
import time

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Cadical153

PORTS = [(l, r) for l in range(4) for r in range(4)]
PIDX = {p: i for i, p in enumerate(PORTS)}
STARS = [[PIDX[(l, r)] for r in range(4)] for l in range(4)] + [[PIDX[(l, r)] for l in range(4)] for r in range(4)]
PMS = [frozenset(PIDX[(l, perm[l])] for l in range(4)) for perm in itertools.permutations(range(4))]
PMSET = set(PMS)


def oa_bijections():
    """All bijections phi: ports -> ports (tuples of length 16) mapping every star to a PM.
    A star maps to a PM iff its four images are pairwise hub-disjoint, so backtrack port by port
    with that pairwise constraint on same-L-star and same-R-star ports."""
    out = []
    phi = [None] * 16
    used = [False] * 16

    def disjoint(a, b):
        return PORTS[a][0] != PORTS[b][0] and PORTS[a][1] != PORTS[b][1]

    def rec(p):
        if p == 16:
            out.append(tuple(phi))
            return
        l, r = PORTS[p]
        for img in range(16):
            if used[img]:
                continue
            ok = True
            for p2 in range(p):
                l2, r2 = PORTS[p2]
                if (l2 == l or r2 == r) and not disjoint(phi[p2], img):
                    ok = False
                    break
            if ok:
                phi[p] = img
                used[img] = True
                rec(p + 1)
                used[img] = False
                phi[p] = None

    rec(0)
    return out


def automorphisms():
    """K4,4 automorphisms as port permutations (1152)."""
    auts = []
    for pl in itertools.permutations(range(4)):
        for pr in itertools.permutations(range(4)):
            for swap in (False, True):
                perm = [0] * 16
                for (l, r), i in PIDX.items():
                    img = (pr[r], pl[l]) if swap else (pl[l], pr[r])
                    perm[i] = PIDX[img]
                auts.append(tuple(perm))
    return auts


def compose_phi(a, phi, b):
    """b o phi o a^-1 as a bijection (a acts on the source, b on the target)."""
    ainv = [0] * 16
    for i, v in enumerate(a):
        ainv[v] = i
    return tuple(b[phi[ainv[x]]] for x in range(16))


def orbit_reps(phis, left_group, right_group):
    phiset = set(phis)
    seen = set()
    reps = []
    for phi in phis:
        if phi in seen:
            continue
        orbit = set()
        for a in left_group:
            for b in right_group:
                orbit.add(compose_phi(a, phi, b))
        assert orbit <= phiset
        seen |= orbit
        reps.append((phi, len(orbit)))
    return reps


class Model120:
    def __init__(self, phi12, phi13):
        self.phi12, self.phi13 = phi12, phi13
        self.inv12 = {v: i for i, v in enumerate(phi12)}
        self.inv13 = {v: i for i, v in enumerate(phi13)}
        self.pool = IDPool()
        self.clauses = []
        v = self.pool.id
        # phi23 permutation variables F[q][r]
        self.F = [[v(("F", q, r)) for r in range(16)] for q in range(16)]
        for q in range(16):
            self.eq1([self.F[q][r] for r in range(16)])
        for r in range(16):
            self.eq1([self.F[q][r] for q in range(16)])
        # OA constraint: image of each star is a PM: for each star S and each hub h of block 3,
        # exactly one q in S with F[q][r] for r in star(h)
        for S in STARS:
            for H in STARS:
                self.eq1([self.F[q][r] for q in S for r in H])
        # edges: E1213[p][p2]: d12[p] ~ d13[p2]; E1223[p][q]: d12[p] ~ d23[q]; E1323[p][q]: d13[p] ~ d23[q]
        self.E1213 = [[v(("E1213", p, p2)) for p2 in range(16)] for p in range(16)]
        self.E1223 = [[v(("E1223", p, q)) for q in range(16)] for p in range(16)]
        self.E1323 = [[v(("E1323", p, q)) for q in range(16)] for p in range(16)]
        for p in range(16):
            self.eqk([self.E1213[p][x] for x in range(16)] + [self.E1223[p][x] for x in range(16)], 2)
            self.eqk([self.E1213[x][p] for x in range(16)] + [self.E1323[p][x] for x in range(16)], 2)
            self.eqk([self.E1223[x][p] for x in range(16)] + [self.E1323[x][p] for x in range(16)], 2)
        self.products = {}
        self.reach_constraints()

    def eq1(self, lits):
        for cl in CardEnc.equals(lits=lits, bound=1, vpool=self.pool, encoding=EncType.seqcounter).clauses:
            self.clauses.append(cl)

    def eqk(self, lits, k):
        for cl in CardEnc.equals(lits=lits, bound=k, vpool=self.pool, encoding=EncType.seqcounter).clauses:
            self.clauses.append(cl)

    def prod(self, a, b):
        key = (min(a, b), max(a, b))
        if key not in self.products:
            w = self.pool.id(("prod",) + key)
            self.clauses += [[-w, a], [-w, b], [w, -a, -b]]
            self.products[key] = w
        return self.products[key]

    def reach_constraints(self):
        phi12, phi13, inv12, inv13 = self.phi12, self.phi13, self.inv12, self.inv13
        F, E1213, E1223, E1323 = self.F, self.E1213, self.E1223, self.E1323
        T = self.pool.id(("true",))
        self.clauses.append([T])
        for star in STARS:              # hub a of block 1
            for q in range(16):         # port of block 2
                lits = []
                for p in star:
                    if phi12[p] == q:
                        lits.append(T)
                    lits.append(E1223[p][q])
                    lits.append(E1213[inv12[q]][p])
                    lits.append(E1323[p][q])
                self.eq1(lits)
            for r in range(16):         # port of block 3
                lits = []
                for p in star:
                    if phi13[p] == r:
                        lits.append(T)
                    lits.append(E1213[p][inv13[r]])
                    for q2 in range(16):
                        lits.append(self.prod(E1323[p][q2], F[q2][r]))
                        lits.append(self.prod(E1223[p][q2], F[q2][r]))
                self.eq1(lits)
        for star in STARS:              # hub b of block 2 (star over block-2 ports q)
            for p in range(16):         # port of block 1
                lits = []
                for q in star:
                    p0 = inv12[q]
                    if p0 == p:
                        lits.append(T)
                    lits.append(E1213[p0][p])
                    lits.append(E1223[p][q])
                    lits.append(E1323[p][q])
                self.eq1(lits)
            for r in range(16):         # port of block 3
                lits = []
                for q in star:
                    p0 = inv12[q]
                    lits.append(F[q][r])
                    lits.append(E1323[inv13[r]][q])
                    lits.append(E1213[p0][inv13[r]])
                    for q2 in range(16):
                        lits.append(self.prod(E1223[p0][q2], F[q2][r]))
                self.eq1(lits)
        for star in STARS:              # hub c of block 3 (star over block-3 ports r)
            for p in range(16):         # port of block 1
                lits = []
                for r in star:
                    p0 = inv13[r]
                    if p0 == p:
                        lits.append(T)
                    lits.append(E1213[p][p0])
                    for q2 in range(16):
                        lits.append(self.prod(F[q2][r], E1223[p][q2]))
                        lits.append(self.prod(F[q2][r], E1323[p][q2]))
                self.eq1(lits)
            for q in range(16):         # port of block 2
                lits = []
                for r in star:
                    p0 = inv13[r]
                    lits.append(F[q][r])
                    for q0 in range(16):
                        lits.append(self.prod(F[q0][r], E1223[inv12[q]][q0]))
                    lits.append(E1213[inv12[q]][p0])
                    lits.append(E1323[p0][q])
                self.eq1(lits)

    def decode(self, model):
        pos = set(l for l in model if l > 0)
        phi23 = [next(r for r in range(16) if self.F[q][r] in pos) for q in range(16)]
        edges = set()
        # vertices: hubs ("H", b, side, i), ports ("P", b, idx), doubles ("D12", p), ("D13", p), ("D23", q)
        for b in range(3):
            for i, (l, r) in enumerate(PORTS):
                edges.add((("P", b, i), ("H", b, "L", l)))
                edges.add((("P", b, i), ("H", b, "R", r)))
        for p in range(16):
            edges.add((("D12", p), ("P", 0, p)))
            edges.add((("D12", p), ("P", 1, self.phi12[p])))
            edges.add((("D13", p), ("P", 0, p)))
            edges.add((("D13", p), ("P", 2, self.phi13[p])))
            edges.add((("D23", p), ("P", 1, p)))
            edges.add((("D23", p), ("P", 2, phi23[p])))
        for p in range(16):
            for x in range(16):
                if self.E1213[p][x] in pos:
                    edges.add((("D12", p), ("D13", x)))
                if self.E1223[p][x] in pos:
                    edges.add((("D12", p), ("D23", x)))
                if self.E1323[p][x] in pos:
                    edges.add((("D13", p), ("D23", x)))
        return phi23, sorted((str(a), str(b)) for a, b in edges)

    def blocking_clause(self, model):
        pos = set(l for l in model if l > 0)
        core = [x for row in self.F for x in row] + [x for M in (self.E1213, self.E1223, self.E1323) for row in M for x in row]
        return [-x for x in core if x in pos]


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
    t0 = time.time()
    phis = oa_bijections()
    auts = automorphisms()
    print("OA bijections (labelled):", len(phis), round(time.time() - t0, 1), "s", file=sys.stderr, flush=True)
    reps12 = orbit_reps(phis, auts, auts)
    print("orbits under Aut x Aut:", [(n) for _, n in reps12], file=sys.stderr, flush=True)
    results = []
    cap = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    for phi12, _ in reps12:
        stab1 = [a for a in auts if any(compose_phi(a, phi12, b) == phi12 for b in auts)]
        # phi13 up to Stab(phi12) on the left and Aut on the right
        reps13 = orbit_reps(phis, stab1, auts)
        print("phi12 rep: |stab1| =", len(stab1), "phi13 reps:", len(reps13), file=sys.stderr, flush=True)
        for phi13, _ in reps13:
            m = Model120(phi12, phi13)
            solver = Cadical153(bootstrap_with=m.clauses)
            t1 = time.time()
            n_sol = 0
            found = []
            while solver.solve():
                model = solver.get_model()
                n_sol += 1
                phi23, edges = m.decode(model)
                chk = bfs_check(edges)
                found.append({"phi23": phi23, "check": chk})
                if chk["diameter"] == 4 and chk["bad_pairs"] == 0:
                    print("!!! DIAMETER-4 GRAPH ON 120 VERTICES", file=sys.stderr, flush=True)
                    with open("runs/graph120_FOUND.json", "w") as f:
                        json.dump({"phi12": phi12, "phi13": phi13, "phi23": phi23, "edges": edges, "check": chk}, f)
                    break
                solver.add_clause(m.blocking_clause(model))
                if n_sol >= cap:
                    break
            rec = {"phi12": phi12, "phi13": phi13, "vars": m.pool.top, "clauses": len(m.clauses),
                   "tight_solutions_found": n_sol, "capped": n_sol >= cap, "solve_s": round(time.time() - t1, 1),
                   "best_bad_pairs": min([f["check"]["bad_pairs"] for f in found], default=None),
                   "examples": found[:3]}
            print(json.dumps({k: rec[k] for k in rec if k not in ("phi12", "phi13", "examples")}), file=sys.stderr, flush=True)
            results.append(rec)
            solver.delete()
    json.dump(results, sys.stdout, indent=1, default=list)


if __name__ == "__main__":
    main()
