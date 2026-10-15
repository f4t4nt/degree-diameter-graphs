"""Structured hunt for n = 110: octahedron + K4,4 + K4,4 blocks with 44 doubles.

Blocks O (octahedron, 6 hubs, 12 ports), A and B (K4,4, 8 hubs, 16 ports each): 66 block vertices,
88 slots. Connectors: 16 doubles d_AB (one per A port, forming an OA bijection A -> B so that the
64 A-B hub pairs are covered exactly once), 12 d_OA and 12 d_OB (every O port holds one of each; the
12 A ports not used by s_A doubles each hold one d_OA, with every A hub carrying exactly three whose
O ports form a perfect matching of the octahedron and every O hub's four d_OA having A ports that
form a perfect matching of K4,4), and two same-block doubles s_A (on a perfect matching of A) and
two s_B. n = 66 + 44 = 110. The hub-to-port budget of this shape has slack, so hub reach is encoded
as at-least-one; solutions are enumerated and checked by BFS for the full diameter condition.

Instance = (phi_AB canonical, choice of the s_A perfect matching and pairing, same for B, up to the
stabiliser symmetries). Variables: psi_A, psi_B (12x12 permutation matrices with the PM constraints)
and the connector edges among the 44 doubles (degree exactly two each).
"""
import collections
import itertools
import json
import sys
import time

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Cadical153

import design120 as D

# octahedron: vertices 0..5, antipodes (0,3), (1,4), (2,5); 12 edges = ports
OCT_EDGES = [(i, j) for i in range(6) for j in range(i + 1, 6) if (i, j) not in ((0, 3), (1, 4), (2, 5))]
OCT_PM = []  # perfect matchings of the octahedron (3 disjoint edges covering 6 vertices)
for trip in itertools.combinations(range(12), 3):
    vs = [v for k in trip for v in OCT_EDGES[k]]
    if len(set(vs)) == 6:
        OCT_PM.append(frozenset(trip))
OCT_STARS = [[k for k, e in enumerate(OCT_EDGES) if v in e] for v in range(6)]
K44_PORTS = D.PORTS  # (l, r)
K44_STARS = D.STARS  # rows then columns, as port indices
K44_PM = D.PMS


class Model110:
    def __init__(self, phi_ab, sA_ports, sA_pairing, sB_ports, sB_pairing):
        """phi_ab: tuple mapping A port index -> B port index. sA_ports: 4 A-port indices (a PM);
        sA_pairing: ((p1,p2),(p3,p4)) the two same-block doubles. Same for B."""
        self.phi = phi_ab
        self.sA, self.sB = sA_pairing, sB_pairing
        self.freeA = [p for p in range(16) if p not in sA_ports]
        self.freeB = [p for p in range(16) if p not in sB_ports]
        self.pool = IDPool()
        self.cl = []
        v = self.pool.id
        self.T = v(("T",))
        self.cl.append([self.T])
        # connectors: 0..15 d_AB (A port p, B port phi[p]); 16..27 d_OA (indexed by O port o);
        # 28..39 d_OB (by O port o); 40,41 s_A; 42,43 s_B
        self.C = list(range(44))
        # psi_A[o][a]: d_OA at O port o holds A port a (a in freeA); psi_B[o][b]
        self.PA = {(o, a): v(("PA", o, a)) for o in range(12) for a in self.freeA}
        self.PB = {(o, b): v(("PB", o, b)) for o in range(12) for b in self.freeB}
        for o in range(12):
            self.eq([self.PA[(o, a)] for a in self.freeA], 1)
            self.eq([self.PB[(o, b)] for b in self.freeB], 1)
        for a in self.freeA:
            self.eq([self.PA[(o, a)] for o in range(12)], 1)
        for b in self.freeB:
            self.eq([self.PB[(o, b)] for o in range(12)], 1)
        # every O hub's four d_OA A-ports form a PM of K4,4: pairwise hub-disjoint
        for star in OCT_STARS:
            for o1, o2 in itertools.combinations(star, 2):
                for a1 in self.freeA:
                    for a2 in self.freeA:
                        if a1 != a2 and not self.k44_disjoint(a1, a2):
                            self.cl.append([-self.PA[(o1, a1)], -self.PA[(o2, a2)]])
                for b1 in self.freeB:
                    for b2 in self.freeB:
                        if b1 != b2 and not self.k44_disjoint(b1, b2):
                            self.cl.append([-self.PB[(o1, b1)], -self.PB[(o2, b2)]])
        # every A hub's d_OA O-ports form a PM of the octahedron: pairwise vertex-disjoint O edges
        for star in K44_STARS:
            fa = [a for a in star if a in self.freeA]
            for a1, a2 in itertools.combinations(fa, 2):
                for o1 in range(12):
                    for o2 in range(12):
                        if o1 != o2 and set(OCT_EDGES[o1]) & set(OCT_EDGES[o2]):
                            self.cl.append([-self.PA[(o1, a1)], -self.PA[(o2, a2)]])
            fb = [b for b in star if b in self.freeB]
            for b1, b2 in itertools.combinations(fb, 2):
                for o1 in range(12):
                    for o2 in range(12):
                        if o1 != o2 and set(OCT_EDGES[o1]) & set(OCT_EDGES[o2]):
                            self.cl.append([-self.PB[(o1, b1)], -self.PB[(o2, b2)]])
        # connector edges
        self.E = {}
        for c in self.C:
            for c2 in self.C:
                if c < c2:
                    self.E[(c, c2)] = self.E[(c2, c)] = v(("E", c, c2))
        for c in self.C:
            self.eq([self.E[(c, c2)] for c2 in self.C if c2 != c], 2)
        self._and = {}
        self.reach()

    @staticmethod
    def k44_disjoint(a1, a2):
        return K44_PORTS[a1][0] != K44_PORTS[a2][0] and K44_PORTS[a1][1] != K44_PORTS[a2][1]

    def eq(self, lits, k):
        for c in CardEnc.equals(lits=lits, bound=k, vpool=self.pool, encoding=EncType.seqcounter).clauses:
            self.cl.append(c)

    def AND(self, *lits):
        lits = tuple(sorted(set(l for l in lits if l != self.T)))
        if not lits:
            return self.T
        if len(lits) == 1:
            return lits[0]
        if lits not in self._and:
            w = self.pool.id(("and",) + lits)
            for l in lits:
                self.cl.append([-w, l])
            self.cl.append([w] + [-l for l in lits])
            self._and[lits] = w
        return self._and[lits]

    # literal: connector c holds A port a / B port b / O port o
    def holdsA(self, c, a):
        if c < 16:
            return self.T if c == a else None
        if 16 <= c < 28:
            return self.PA.get((c - 16, a))
        if c in (40, 41):
            return self.T if a in self.sA[c - 40] else None
        return None

    def holdsB(self, c, b):
        if c < 16:
            return self.T if self.phi[c] == b else None
        if 28 <= c < 40:
            return self.PB.get((c - 28, b))
        if c in (42, 43):
            return self.T if b in self.sB[c - 42] else None
        return None

    def holdsO(self, c, o):
        if 16 <= c < 28:
            return self.T if c - 16 == o else None
        if 28 <= c < 40:
            return self.T if c - 28 == o else None
        return None

    def reach(self):
        """Hub-to-port at-least-one for all cross-block (hub, port) pairs: a-p-c-q or a-p-c-c'-q."""
        holds = {"A": self.holdsA, "B": self.holdsB, "O": self.holdsO}
        ports_of_hub = {"A": [(a, "A") for a in range(16)], "B": [(b, "B") for b in range(16)], "O": [(o, "O") for o in range(12)]}
        hubs = {"A": K44_STARS, "B": K44_STARS, "O": OCT_STARS}
        nports = {"A": 16, "B": 16, "O": 12}
        self.n_reach = 0
        for blk, stars in hubs.items():
            for star in stars:
                for tgt in ("A", "B", "O"):
                    if tgt == blk:
                        continue
                    for q in range(nports[tgt]):
                        lits = []
                        for p in star:
                            for c in self.C:
                                hp = holds[blk](c, p)
                                if hp is None:
                                    continue
                                hq = holds[tgt](c, q)
                                if hq is not None:
                                    lits.append(self.AND(hp, hq))
                                for c2 in self.C:
                                    if c2 == c:
                                        continue
                                    hq2 = holds[tgt](c2, q)
                                    if hq2 is not None:
                                        lits.append(self.AND(hp, self.E[(c, c2)], hq2))
                        lits = [l for l in lits if l is not None]
                        if self.T in lits:
                            continue
                        self.cl.append(lits)
                        self.n_reach += 1

    def decode(self, model):
        pos = set(l for l in model if l > 0)
        edges = set()
        for a, (l, r) in enumerate(K44_PORTS):
            edges.add((("P", "A", a), ("H", "A", "L", l))); edges.add((("P", "A", a), ("H", "A", "R", r)))
            edges.add((("P", "B", a), ("H", "B", "L", l))); edges.add((("P", "B", a), ("H", "B", "R", r)))
        for o, (u, w) in enumerate(OCT_EDGES):
            edges.add((("P", "O", o), ("H", "O", u))); edges.add((("P", "O", o), ("H", "O", w)))
        for c in self.C:
            for a in range(16):
                h = self.holdsA(c, a)
                if h is not None and (h == self.T or h in pos):
                    edges.add((("C", c), ("P", "A", a)))
            for b in range(16):
                h = self.holdsB(c, b)
                if h is not None and (h == self.T or h in pos):
                    edges.add((("C", c), ("P", "B", b)))
            for o in range(12):
                h = self.holdsO(c, o)
                if h is not None and (h == self.T or h in pos):
                    edges.add((("C", c), ("P", "O", o)))
            for c2 in self.C:
                if c < c2 and self.E[(c, c2)] in pos:
                    edges.add((("C", c), ("C", c2)))
        return sorted((str(x), str(y)) for x, y in edges)

    def blocking(self, model):
        pos = set(l for l in model if l > 0)
        core = list(self.PA.values()) + list(self.PB.values()) + [v for (c, c2), v in self.E.items() if c < c2]
        return [-x for x in core if x in pos]


def main():
    cap = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    phis = D.oa_bijections()
    auts = D.automorphisms()
    phi = D.orbit_reps(phis, auts, auts)[0][0]
    # s_A perfect matchings up to the left stabiliser of phi; pairings of a PM into two pairs: 3
    stab_left = [a for a in auts if any(D.compose_phi(a, phi, b) == phi for b in auts)]
    pm_reps = []
    seen = set()
    for pm in K44_PM:
        if pm in seen:
            continue
        orbit = {frozenset(a[p] for p in pm) for a in stab_left}
        seen |= orbit
        pm_reps.append(pm)
    print("phi stabiliser (left):", len(stab_left), "s_A PM representatives:", len(pm_reps), file=sys.stderr, flush=True)
    results = []
    t0 = time.time()
    for pmA in pm_reps:
        pa = sorted(pmA)
        for pairA in [((pa[0], pa[1]), (pa[2], pa[3])), ((pa[0], pa[2]), (pa[1], pa[3])), ((pa[0], pa[3]), (pa[1], pa[2]))]:
            for pmB in K44_PM:
                pb = sorted(pmB)
                for pairB in [((pb[0], pb[1]), (pb[2], pb[3])), ((pb[0], pb[2]), (pb[1], pb[3])), ((pb[0], pb[3]), (pb[1], pb[2]))]:
                    m = Model110(phi, pa, pairA, pb, pairB)
                    s = Cadical153(bootstrap_with=m.cl)
                    n_sol = 0
                    best = None
                    t1 = time.time()
                    while s.solve():
                        n_sol += 1
                        model = s.get_model()
                        edges = m.decode(model)
                        chk = D.bfs_check(edges)
                        if best is None or chk["bad_pairs"] < best["bad_pairs"]:
                            best = chk
                        if chk["bad_pairs"] == 0 and chk["diameter"] == 4:
                            print("!!! VERIFIED DIAMETER-4 GRAPH n =", chk["n"], file=sys.stderr, flush=True)
                            json.dump({"edges": edges, "check": chk, "pmA": pa, "pairA": pairA, "pmB": pb, "pairB": pairB},
                                      open("runs/graph110_FOUND.json", "w"))
                            break
                        s.add_clause(m.blocking(model))
                        if n_sol >= cap:
                            break
                    s.delete()
                    rec = {"pmA": pa, "pairA": pairA, "pmB": pb, "pairB": pairB, "vars": m.pool.top, "clauses": len(m.cl),
                           "reach_clauses": m.n_reach, "hub_reach_solutions": n_sol, "capped": n_sol >= cap,
                           "best_check": best, "solve_s": round(time.time() - t1, 1)}
                    results.append(rec)
                    print(json.dumps({k: rec[k] for k in ("pmA", "pairA", "pmB", "pairB", "hub_reach_solutions", "capped", "best_check", "solve_s")}),
                          file=sys.stderr, flush=True)
                    if best is not None and best["bad_pairs"] == 0:
                        json.dump(results, sys.stdout, indent=1)
                        return
    print("total", round(time.time() - t0, 1), "s", file=sys.stderr)
    json.dump(results, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
