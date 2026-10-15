"""Typed SAT model for three subdivided K4,4 blocks plus connectors, with slack (no exactness).

Input shape (shapes106.py): u units and a list of gadgets. A unit is two adjacent triples with
exactly two ports per block in total (types free: pure (1,1,1)+(1,1,1) or bent complementary pairs).
A gadget is a typed connector component: connector profiles x(c) and connector edges fixed.

Variables A[c][q] (connector c holds port q). Cardinalities: a gadget connector holds exactly x_b(c)
ports of block b; a unit member holds 3 ports and the unit 2 per block; every port has exactly two
holders. Derived (both directions where cheap): U[c][a] (c holds a port of hub a), Z[c][q] (c within
2 of q: own port or a connector neighbour's), S[c][c'] (c and c' share a port), R3[c][q], R4[c][q]
(c within 3 / 4 of q), with every walk shape of TEMPLATE_CALCULUS.md section 3 (ports are adjacent
only to hubs and connectors, hubs only to ports, no port-port edges, hubs saturated):

  hub a - foreign port q: a-p-c-q, a-p-c-c'-q               -> OR_c U[c][a] & Z[c][q]
  hub a - hub b (other block): a-p-c-q-b                      -> OR_c U[c][a] & U[c][b]
  hub a - connector c without a block(a) port: a-p-c0-c, a-p-c0-c1-c, a-p-c0-q-c
  port p - port q (other block): p-c-(<=3 more)-q, p-a-p'-c-q -> OR_c A[c][p] & R3[c][q] | HP3
  port p - connector c: R4[c][p]
  connector c - c' at connector-graph distance > 4: via shared ports (S) and hubs (U & U)
  within a block everything is within 4 (hub-hub 2/4, hub-port <= 3, port-port 2/4), and a connector
  holding a block-i port is within 4 of every hub of block i (a-p-b-p'-c).

Redundant: port lemma inequality (section 7 in inequality form): for every port q of block j and
pair of potential holders with f_i(h1) + f_i(h2) < 4 for some i != j, not both (f(c) = carried
vector, fixed by the shape; unit members have f = (2,2,2)); for two units sharing a port, their
four block-i ports (i != j) carry exactly one route to every hub of block i (4 ports, 8 hub
incidences, 8 hubs), encoded per hub. Symmetry breaking: per block the first gadget connector
with a port there (else unit 0) holds port 0 of the block, a second port of that connector is
(0,1) or (1,1); interchangeable units are ordered by their smallest block-0 port; inside a unit
the member holding the unit's smallest block-0 port is member 0.

Usage:
  typed_slack_sat.py run SHAPES.json [--idx=i,j,...] [--cap=SECONDS] [--out=runs/NAME.json]
                     [--full] [--notight] [--noexact] [--nolemma]
      default: stage 1 hub-to-port/hub-to-hub model, stage 2 full model on SAT; --full: full model
      directly (used for the 106 batch: the hub model alone does not refute the triangle shapes);
      tight mode is switched on automatically for budget-768 shapes unless --notight.
  typed_slack_sat.py validate FILE.edges          known 104: own assignment SAT, free model SAT + BFS
  typed_slack_sat.py control [hub|full|tight|tight_hub|noexact|nounits|nolemma] [--cap=S]
                                                  105 cross triangle (must be UNSAT, sat105 agrees)

Caveat: python-sat 1.9.dev15 does not support Cadical153.interrupt(), so --cap is not enforced per
instance (the timer thread raises NotImplementedError and the solve runs to completion); bound runs
with an outer `timeout`, or switch to conf_budget.
"""
import collections
import hashlib
import itertools
import json
import pathlib
import sys
import threading
import time

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Cadical153

HERE = pathlib.Path(__file__).resolve().parent
NB = 3
PORTS = [(b, l, r) for b in range(NB) for l in range(4) for r in range(4)]
PIDX = {p: i for i, p in enumerate(PORTS)}
HUBS = [(b, s, k) for b in range(NB) for s in range(2) for k in range(4)]  # s=0: L side, s=1: R side
HIDX = {h: i for i, h in enumerate(HUBS)}


def ports_of_hub(h):
    b, s, k = h
    return [PIDX[(b, k, r)] for r in range(4)] if s == 0 else [PIDX[(b, l, k)] for l in range(4)]


def hubs_of_port(pi):
    b, l, r = PORTS[pi]
    return [HIDX[(b, 0, l)], HIDX[(b, 1, r)]]


PORTS_OF_BLOCK = {b: [i for i, p in enumerate(PORTS) if p[0] == b] for b in range(NB)}
HUBS_OF_BLOCK = {b: [i for i, h in enumerate(HUBS) if h[0] == b] for b in range(NB)}
POH = [ports_of_hub(h) for h in HUBS]
HOP = [hubs_of_port(i) for i in range(len(PORTS))]


class Model:
    def __init__(self, units, gadgets, exact_units=True, lemma=True, symbreak=True, full=True, exact_all=True, tight=False):
        """full=False: hub-to-port and hub-to-hub constraints only (plus cardinalities, redundant and
        symmetry clauses); a model that is UNSAT without the port/connector distance constraints is
        UNSAT with them, and a SAT answer is then re-run with full=True.
        tight=True (budget exactly 768 shapes only): every hub reaches every foreign port exactly once,
        so a connector holds at most one port of a hub and the route literals of each (hub, foreign
        port) pair are at most one."""
        self.tight = tight
        self.pool = IDPool()
        self.cl = []
        v = self.pool.id
        # connectors: gadget connectors first, then unit members
        self.types = []      # fixed profile or None (unit member)
        self.unit_of = []    # unit index or None
        self.edges = set()
        self.gadget_of = []
        for gi, g in enumerate(gadgets):
            base = len(self.types)
            for t in g["types"]:
                self.types.append(tuple(t))
                self.unit_of.append(None)
                self.gadget_of.append(gi)
            for a, b in g["edges"]:
                self.edges.add((base + min(a, b), base + max(a, b)))
        self.units = []
        for ui in range(units):
            h1, h2 = len(self.types), len(self.types) + 1
            self.types += [None, None]
            self.unit_of += [ui, ui]
            self.gadget_of += [None, None]
            self.edges.add((h1, h2))
            self.units.append((h1, h2))
        NC = self.nc = len(self.types)
        self.nbr = {c: set() for c in range(NC)}
        for a, b in self.edges:
            self.nbr[a].add(b)
            self.nbr[b].add(a)
        # carried vectors f (fixed for every connector)
        self.f = []
        for c in range(NC):
            if self.types[c] is None:
                self.f.append((2, 2, 2))
            else:
                self.f.append(tuple(self.types[c][j] + sum(self.types[d][j] for d in self.nbr[c]) for j in range(3)))
        # possible ports per connector
        self.poss = []
        for c in range(NC):
            if self.types[c] is None:
                self.poss.append(set(range(len(PORTS))))
            else:
                self.poss.append({i for i in range(len(PORTS)) if self.types[c][PORTS[i][0]] > 0})
        self.A = {}
        for c in range(NC):
            for i in self.poss[c]:
                self.A[(c, i)] = v(("A", c, i))
        self.holders = {i: [c for c in range(NC) if i in self.poss[c]] for i in range(len(PORTS))}
        # cardinalities
        for c in range(NC):
            if self.types[c] is None:
                self.eq([self.A[(c, i)] for i in range(len(PORTS))], 3)
            else:
                for b in range(NB):
                    if self.types[c][b]:
                        self.eq([self.A[(c, i)] for i in PORTS_OF_BLOCK[b]], self.types[c][b])
        for (h1, h2) in self.units:
            for b in range(NB):
                self.eq([self.A[(h1, i)] for i in PORTS_OF_BLOCK[b]] + [self.A[(h2, i)] for i in PORTS_OF_BLOCK[b]], 2)
        for i in range(len(PORTS)):
            self.eq([self.A[(c, i)] for c in self.holders[i]], 2)
        # U[c][a]
        self.U = {}
        for c in range(NC):
            for a in range(len(HUBS)):
                lits = [self.A[(c, i)] for i in POH[a] if (c, i) in self.A]
                if not lits:
                    continue
                u = v(("U", c, a))
                self.U[(c, a)] = u
                self.cl.append([-u] + lits)
                for l in lits:
                    self.cl.append([u, -l])
                if tight:
                    for l1, l2 in itertools.combinations(lits, 2):
                        self.cl.append([-l1, -l2])
        # Z[c][q]: own port or a neighbour's
        self.Z = {}
        for c in range(NC):
            for i in range(len(PORTS)):
                lits = [self.A[(d, i)] for d in sorted(self.nbr[c] | {c}) if (d, i) in self.A]
                if not lits:
                    continue
                z = v(("Z", c, i))
                self.Z[(c, i)] = z
                self.cl.append([-z] + lits)
                for l in lits:
                    self.cl.append([z, -l])
        if full:
            self.distance_constraints()
        # hub - foreign port, and HP3[a][q] (a reaches q within 3: a-p-c-q)
        self.HP3 = {}
        for a in range(len(HUBS)):
            for i in range(len(PORTS)):
                if PORTS[i][0] == HUBS[a][0]:
                    continue
                terms, terms3 = [], []
                for c in range(NC):
                    if (c, a) not in self.U:
                        continue
                    if (c, i) in self.Z:
                        w = v(("W", c, a, i))
                        self.cl += [[-w, self.U[(c, a)]], [-w, self.Z[(c, i)]]]
                        terms.append(w)
                    if (c, i) in self.A:
                        w3 = v(("W3", c, a, i))
                        self.cl += [[-w3, self.U[(c, a)]], [-w3, self.A[(c, i)]]]
                        terms3.append(w3)
                self.cl.append(terms)
                h3 = v(("HP3", a, i))
                self.HP3[(a, i)] = h3
                self.cl.append([-h3] + terms3)
                if tight:
                    # exactly one route a-p-c-q / a-p-c-c''-q: one literal per (slot connector, carrier)
                    routes = []
                    for c in range(NC):
                        if (c, a) not in self.U:
                            continue
                        for c2 in sorted(self.nbr[c] | {c}):
                            if (c2, i) in self.A:
                                r = v(("RT", c, c2, a, i))
                                self.cl += [[-r, self.U[(c, a)]], [-r, self.A[(c2, i)]], [r, -self.U[(c, a)], -self.A[(c2, i)]]]
                                routes.append(r)
                    for cl in CardEnc.atmost(lits=routes, bound=1, vpool=self.pool, encoding=EncType.seqcounter).clauses:
                        self.cl.append(cl)
        # hub - hub across blocks
        for a in range(len(HUBS)):
            for b in range(a + 1, len(HUBS)):
                if HUBS[a][0] == HUBS[b][0]:
                    continue
                terms = []
                for c in range(NC):
                    if (c, a) in self.U and (c, b) in self.U:
                        w = v(("HH", c, a, b))
                        self.cl += [[-w, self.U[(c, a)]], [-w, self.U[(c, b)]]]
                        terms.append(w)
                self.cl.append(terms)
        if full:
            self.port_connector_constraints()
        self.n_core = len(self.cl)
        if lemma:
            # port lemma, inequality form
            for i in range(len(PORTS)):
                j = PORTS[i][0]
                for c, d in itertools.combinations(self.holders[i], 2):
                    if any(self.f[c][k] + self.f[d][k] < 4 for k in range(3) if k != j):
                        self.cl.append([-self.A[(c, i)], -self.A[(d, i)]])
        if exact_all:
            self.exactness_gadget_pairs()
        if exact_units and len(self.units) >= 2:
            # two units sharing a port of block j carry exactly one route to each hub of block i != j;
            # the two members of one unit never share a port (their 8 routes reach at most 4 hubs)
            UP, UH1, UH2 = {}, {}, {}
            for ui, (h1, h2) in enumerate(self.units):
                for i in range(len(PORTS)):
                    p = v(("UP", ui, i))
                    UP[(ui, i)] = p
                    self.cl += [[-p, self.A[(h1, i)], self.A[(h2, i)]], [p, -self.A[(h1, i)]], [p, -self.A[(h2, i)]]]
                    self.cl.append([-self.A[(h1, i)], -self.A[(h2, i)]])
                for a in range(len(HUBS)):
                    lits = [UP[(ui, i)] for i in POH[a]]
                    g1 = v(("UH1", ui, a))
                    g2 = v(("UH2", ui, a))
                    UH1[(ui, a)], UH2[(ui, a)] = g1, g2
                    self.cl.append([-g1] + lits)
                    for l in lits:
                        self.cl.append([g1, -l])
                    # g2 <-> at least two of lits (a unit has at most two ports per block, so a hub
                    # with two of them means both block ports share this hub)
                    for l1, l2 in itertools.combinations(lits, 2):
                        self.cl.append([g2, -l1, -l2])
                    self.cl.append([-g2] + [v(("UH2p", ui, a, k)) for k in range(len(lits) * (len(lits) - 1) // 2)])
                    for k, (l1, l2) in enumerate(itertools.combinations(lits, 2)):
                        self.cl += [[-v(("UH2p", ui, a, k)), l1], [-v(("UH2p", ui, a, k)), l2]]
            for u1, u2 in itertools.combinations(range(len(self.units)), 2):
                for i in range(len(PORTS)):
                    j = PORTS[i][0]
                    sh = [-UP[(u1, i)], -UP[(u2, i)]]
                    for a in range(len(HUBS)):
                        if HUBS[a][0] == j:
                            continue
                        self.cl.append(sh + [-UH1[(u1, a)], -UH1[(u2, a)]])
                        self.cl.append(sh + [-UH2[(u1, a)]])
                        self.cl.append(sh + [-UH2[(u2, a)]])
                        self.cl.append(sh + [UH1[(u1, a)], UH1[(u2, a)]])
        self.n_redundant = len(self.cl) - self.n_core
        if symbreak:
            self.symmetry_breaking()
        self.n_sym = len(self.cl) - self.n_core - self.n_redundant

    def distance_constraints(self):
        """S[c][c'] (shared port), R3[c][q], R4[c][q] (c within 3 / 4 of q); R4 must hold everywhere."""
        v = self.pool.id
        NC = self.nc
        self.S = {}
        for c in range(NC):
            for d in range(c + 1, NC):
                common = sorted(self.poss[c] & self.poss[d])
                if not common:
                    continue
                s = v(("S", c, d))
                self.S[(c, d)] = self.S[(d, c)] = s
                terms = []
                for i in common:
                    t = v(("SP", c, d, i))
                    self.cl += [[-t, self.A[(c, i)]], [-t, self.A[(d, i)]], [t, -self.A[(c, i)], -self.A[(d, i)]]]
                    terms.append(t)
                self.cl.append([-s] + terms)
                for t in terms:
                    self.cl.append([s, -t])
        # R3[c][q] (soundness direction only; R3 may be false even if c is within 3)
        self.R3 = {}
        for c in range(NC):
            for i in range(len(PORTS)):
                terms = []
                if (c, i) in self.Z:
                    terms.append(self.Z[(c, i)])
                for d in self.nbr[c]:
                    if (d, i) in self.Z:
                        terms.append(self.Z[(d, i)])
                for a in HOP[i]:
                    if (c, a) in self.U:
                        terms.append(self.U[(c, a)])
                for d in range(NC):
                    if d != c and (c, d) in self.S and (d, i) in self.A:
                        y = v(("Y3", c, d, i))
                        self.cl += [[-y, self.S[(c, d)]], [-y, self.A[(d, i)]]]
                        terms.append(y)
                r = v(("R3", c, i))
                self.R3[(c, i)] = r
                self.cl.append([-r] + terms)
        # R4[c][q] must hold for every connector and port
        for c in range(NC):
            for i in range(len(PORTS)):
                terms = [self.R3[(c, i)]] + [self.R3[(d, i)] for d in self.nbr[c]]
                for d in range(NC):
                    if d != c and (c, d) in self.S and (d, i) in self.Z:
                        y = v(("Y4", c, d, i))
                        self.cl += [[-y, self.S[(c, d)]], [-y, self.Z[(d, i)]]]
                        terms.append(y)
                self.cl.append(terms)

    def port_connector_constraints(self):
        """hub - connector (connector without a block(a) port), port - port across blocks,
        connector - connector at connector-graph distance > 4."""
        v = self.pool.id
        NC = self.nc
        n2 = {c: set(itertools.chain.from_iterable(self.nbr[d] for d in self.nbr[c])) | self.nbr[c] for c in range(NC)}
        for a in range(len(HUBS)):
            ba = HUBS[a][0]
            for c in range(NC):
                if self.types[c] is None or self.types[c][ba] > 0:
                    continue  # a connector with a block(a) port is within 4 of a (a-p-b-p'-c)
                terms = [self.U[(d, a)] for d in n2[c] if d != c and (d, a) in self.U]
                for d in range(NC):
                    if d != c and (d, a) in self.U and (c, d) in self.S:
                        w = v(("HC", c, d, a))
                        self.cl += [[-w, self.U[(d, a)]], [-w, self.S[(c, d)]]]
                        terms.append(w)
                self.cl.append(terms)
        for i in range(len(PORTS)):
            for j in range(len(PORTS)):
                if PORTS[i][0] >= PORTS[j][0]:
                    continue
                terms = [self.HP3[(a, j)] for a in HOP[i]]
                for c in self.holders[i]:
                    w = v(("PP", c, i, j))
                    self.cl += [[-w, self.A[(c, i)]], [-w, self.R3[(c, j)]]]
                    terms.append(w)
                self.cl.append(terms)
        dist = self.connector_distances()
        for c in range(NC):
            for d in range(c + 1, NC):
                if dist[c].get(d, 99) <= 4:
                    continue
                terms = []
                if (c, d) in self.S:
                    terms.append(self.S[(c, d)])
                for c1 in n2[c]:
                    if (c1, d) in self.S:
                        terms.append(self.S[(c1, d)])
                for d1 in n2[d]:
                    if (c, d1) in self.S:
                        terms.append(self.S[(c, d1)])
                for c1 in self.nbr[c]:
                    for d1 in self.nbr[d]:
                        if (c1, d1) in self.S:
                            terms.append(self.S[(c1, d1)])
                for a in range(len(HUBS)):
                    if (c, a) in self.U and (d, a) in self.U:
                        w = v(("CCU", c, d, a))
                        self.cl += [[-w, self.U[(c, a)]], [-w, self.U[(d, a)]]]
                        terms.append(w)
                for e in range(NC):
                    if e != c and e != d and (c, e) in self.S and (e, d) in self.S:
                        w = v(("CCS", c, d, e))
                        self.cl += [[-w, self.S[(c, e)]], [-w, self.S[(e, d)]]]
                        terms.append(w)
                self.cl.append(terms)

    def exactness_gadget_pairs(self):
        """Exactly one route per hub at tight holder pairs involving a gadget connector.

        At a port q of block j held by h1, h2 the routes from the hubs of block i != j are
        a-p-c-q / a-p-c-c'-q with c in the multiset N[h1] + N[h2]; their number is
        2 (f_i(h1) + f_i(h2)). When f_i(h1) + f_i(h2) = 4 (a constant of the shape for every pair of
        connectors) the 8 routes must cover the 8 hubs, so no hub is reached twice: no two members
        of the multiset carry a port at the same hub, a member listed twice carries no block-i port,
        and no member holds two ports of one hub. NH[h][a] = some c in N[h] holds a port of hub a,
        NH2[h][a] = N[h] carries two incidences at a (two members, or one member with two ports);
        both in the forcing direction. Unit-unit pairs are handled at unit level elsewhere."""
        v = self.pool.id
        NC = self.nc
        gad = [c for c in range(NC) if self.types[c] is not None]
        if not gad:
            return
        # HI2[c][a]: c holds two ports of hub a (gadget connectors with >= 2 ports in that block, unit members)
        HI2 = {}
        for c in range(NC):
            for a in range(len(HUBS)):
                lits = [self.A[(c, i)] for i in POH[a] if (c, i) in self.A]
                if len(lits) < 2 or (self.types[c] is not None and self.types[c][HUBS[a][0]] < 2):
                    continue
                g2 = v(("HI2", c, a))
                HI2[(c, a)] = g2
                for l1, l2 in itertools.combinations(lits, 2):
                    self.cl.append([g2, -l1, -l2])
        NH, NH2 = {}, {}
        for h in range(NC):
            closed = sorted(self.nbr[h] | {h})
            for a in range(len(HUBS)):
                us = [self.U[(c, a)] for c in closed if (c, a) in self.U]
                if not us:
                    continue
                n1 = v(("NH", h, a))
                NH[(h, a)] = n1
                self.cl.append([-n1] + us)
                for l in us:
                    self.cl.append([n1, -l])
                n2 = v(("NH2", h, a))
                NH2[(h, a)] = n2
                for l1, l2 in itertools.combinations(us, 2):
                    self.cl.append([n2, -l1, -l2])
                for c in closed:
                    if (c, a) in HI2:
                        self.cl.append([n2, -HI2[(c, a)]])
        for i in range(len(PORTS)):
            j = PORTS[i][0]
            for h1, h2 in itertools.combinations(self.holders[i], 2):
                if self.types[h1] is None and self.types[h2] is None:
                    continue
                for b in range(NB):
                    if b == j or self.f[h1][b] + self.f[h2][b] != 4:
                        continue
                    both = (self.nbr[h1] | {h1}) & (self.nbr[h2] | {h2})
                    pre = [-self.A[(h1, i)], -self.A[(h2, i)]]
                    for a in HUBS_OF_BLOCK[b]:
                        if (h1, a) in NH and (h2, a) in NH:
                            self.cl.append(pre + [-NH[(h1, a)], -NH[(h2, a)]])
                        for h in (h1, h2):
                            if (h, a) in NH2:
                                self.cl.append(pre + [-NH2[(h, a)]])
                        for c in both:
                            if (c, a) in self.U:
                                self.cl.append(pre + [-self.U[(c, a)]])

    def connector_distances(self):
        dist = {}
        for c in range(self.nc):
            d = {c: 0}
            fr = [c]
            while fr:
                nx_ = []
                for x in fr:
                    for y in self.nbr[x]:
                        if y not in d:
                            d[y] = d[x] + 1
                            nx_.append(y)
                fr = nx_
            dist[c] = d
        return dist

    def symmetry_breaking(self):
        v = self.pool.id
        anchored_unit0 = False
        for b in range(NB):
            p0 = PIDX[(b, 0, 0)]
            gad = [c for c in range(self.nc) if self.types[c] is not None and self.types[c][b] > 0]
            if gad:
                c0 = gad[0]
                self.cl.append([self.A[(c0, p0)]])
                if self.types[c0][b] >= 2:
                    # stabiliser of (0,0): S3 x S3 on the other hubs and the side swap; a second port
                    # either shares a hub with (0,0) -> (0,1), or is disjoint -> (1,1)
                    self.cl.append([self.A[(c0, PIDX[(b, 0, 1)])], self.A[(c0, PIDX[(b, 1, 1)])]])
            elif self.units:
                h1, h2 = self.units[0]
                self.cl.append([self.A[(h1, p0)], self.A[(h2, p0)]])
                anchored_unit0 = True
        # interchangeable units ordered by smallest block-0 port; unit 0 excluded when it is an anchor
        order = list(range(1 if anchored_unit0 else 0, len(self.units)))
        if len(order) >= 2:
            F = {}
            for ui in order:
                h1, h2 = self.units[ui]
                for k, i in enumerate(PORTS_OF_BLOCK[0]):
                    f = v(("F", ui, k))
                    F[(ui, k)] = f
                    lits = [self.A[(h1, i2)] for i2 in PORTS_OF_BLOCK[0][:k + 1]] + [self.A[(h2, i2)] for i2 in PORTS_OF_BLOCK[0][:k + 1]]
                    self.cl.append([-f] + lits)
                    for l in lits:
                        self.cl.append([f, -l])
            for ua, ub in zip(order, order[1:]):
                for k in range(16):
                    self.cl.append([-F[(ub, k)], F[(ua, k)]])
        # inside a unit: member 0 holds the unit's smallest block-0 port
        for (h1, h2) in self.units:
            for k, i in enumerate(PORTS_OF_BLOCK[0]):
                self.cl.append([-self.A[(h2, i)]] + [self.A[(h1, i2)] for i2 in PORTS_OF_BLOCK[0][:k]] + [self.A[(h2, i2)] for i2 in PORTS_OF_BLOCK[0][:k]])

    def eq(self, lits, k):
        for cl in CardEnc.equals(lits=lits, bound=k, vpool=self.pool, encoding=EncType.seqcounter).clauses:
            self.cl.append(cl)

    def decode(self, model):
        """Integer-labelled edge list: hubs 0..23, ports 24..71, connectors 72.."""
        pos = set(l for l in model if l > 0)
        edges = set()
        for i in range(len(PORTS)):
            for a in HOP[i]:
                edges.add((a, 24 + i))
        for (c, i), lit in self.A.items():
            if lit in pos:
                edges.add((24 + i, 72 + c))
        for a, b in self.edges:
            edges.add((72 + a, 72 + b))
        return sorted(edges)


def bfs_check(edges):
    adj = collections.defaultdict(set)
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    verts = sorted(adj)
    bad = 0
    diam = 0
    for s in verts:
        dist = {s: 0}
        fr = [s]
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
        bad += sum(1 for u, d in dist.items() if d > 4 and u > s)
    return {"n": len(verts), "degrees": dict(collections.Counter(len(adj[x]) for x in verts)), "diameter": diam, "bad_pairs": bad}


def solve_capped(clauses, cap, assumptions=()):
    s = Cadical153(bootstrap_with=clauses)
    timer = threading.Timer(cap, s.interrupt) if cap else None
    if timer:
        timer.start()
    t0 = time.time()
    r = s.solve_limited(assumptions=list(assumptions), expect_interrupt=True)
    el = time.time() - t0
    if timer:
        timer.cancel()
    model = s.get_model() if r else None
    s.delete()
    return r, el, model


def run_shape(shape, idx, cap, exact_units, lemma, outdir, full=False, allow_tight=True, exact_all=True):
    """Stage 1: hub-to-port / hub-to-hub model (full=False). UNSAT there is UNSAT for the graph
    problem. SAT: stage 2 with the full distance model under the remaining cap; a stage-2 model is
    decoded and BFS-checked independently."""
    t0 = time.time()
    tight = shape.get("budget") == 768 and allow_tight
    m = Model(shape["units"], shape["gadgets"], exact_units=exact_units, lemma=lemma, full=full, tight=tight, exact_all=exact_all)
    rec = {"idx": idx, "units": shape["units"], "potential": shape["potential"], "budget": shape.get("budget"),
           "gadgets": [(g["types"], g["edges"]) for g in shape["gadgets"]], "tight": tight,
           "vars": m.pool.top, "clauses": len(m.cl), "core": m.n_core, "redundant": m.n_redundant, "sym": m.n_sym,
           "build_s": round(time.time() - t0, 1), "stage1_full": full}
    r, el, model = solve_capped(m.cl, cap)
    rec["stage1"] = {True: "SAT", False: "UNSAT", None: "TIMEOUT"}[r]
    rec["stage1_s"] = round(el, 1)
    rec["result"] = rec["stage1"]
    if r and not full:
        rec["stage1_check"] = bfs_check(m.decode(model))
        t1 = time.time()
        m2 = Model(shape["units"], shape["gadgets"], exact_units=exact_units, lemma=lemma, full=True, tight=tight)
        rec["stage2_vars"], rec["stage2_clauses"], rec["stage2_build_s"] = m2.pool.top, len(m2.cl), round(time.time() - t1, 1)
        r, el, model = solve_capped(m2.cl, max(60, cap - el))
        rec["stage2"] = {True: "SAT", False: "UNSAT", None: "TIMEOUT"}[r]
        rec["stage2_s"] = round(el, 1)
        rec["result"] = rec["stage2"]
        m = m2
    if r:
        edges = m.decode(model)
        rec["check"] = bfs_check(edges)
        n = rec["check"]["n"]
        path = outdir / f"solution_{n}_idx{idx}.edges"
        path.write_text("".join(f"{a} {b}\n" for a, b in edges))
        rec["edges_file"] = str(path)
        rec["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return rec


def main_run(args):
    shapes = json.load(open(args[0]))["shapes"]
    idxs = next((a.split("=")[1] for a in args if a.startswith("--idx=")), None)
    idxs = [int(x) for x in idxs.split(",")] if idxs else list(range(len(shapes)))
    cap = float(next((a.split("=")[1] for a in args if a.startswith("--cap=")), 600))
    out = next((a.split("=")[1] for a in args if a.startswith("--out=")), None)
    exact_units = "--noexact" not in args
    lemma = "--nolemma" not in args
    full = "--full" in args
    allow_tight = "--notight" not in args
    exact_all = "--noexact" not in args
    outdir = HERE / "runs"
    results = []
    for i in idxs:
        rec = run_shape(shapes[i], i, cap, exact_units, lemma, outdir, full=full, allow_tight=allow_tight, exact_all=exact_all)
        print(json.dumps(rec), file=sys.stderr, flush=True)
        results.append(rec)
        if out:
            json.dump(results, open(out, "w"), indent=1)
    if not out:
        json.dump(results, sys.stdout, indent=1)


def validate(path):
    sys.path.insert(0, str(HERE.parent))
    import verify

    edges = verify.load_edges_file(path)
    verts, adj, problems = verify.build(edges)
    assert not problems
    blocks = verify.k44_blocks(verts, adj)
    decs = verify.decompositions(blocks, len(verts))
    assert len(decs) == 1
    trip = decs[0]
    hub_of, port_of = {}, {}
    for b, key in enumerate(trip):
        L, R, pmap = blocks[key]
        for k, l in enumerate(L):
            hub_of[l] = HIDX[(b, 0, k)]
        for k, r in enumerate(R):
            hub_of[r] = HIDX[(b, 1, k)]
        for (l, r), q in pmap.items():
            port_of[q] = PIDX[(b, L.index(l), R.index(r))]
    conns = [x for x in verts if x not in hub_of and x not in port_of]
    cadj = {c: [u for u in adj[c] if u in set(conns)] for c in conns}
    assert all(len(cadj[c]) == 1 for c in conns)
    pairs, seen = [], set()
    for c in conns:
        if c not in seen:
            pairs.append((c, cadj[c][0]))
            seen |= {c, cadj[c][0]}
    rec = {"file": path, "units": len(pairs)}
    # (1) free-type units, no symmetry breaking, own assignment under assumptions
    m = Model(len(pairs), [], symbreak=False)
    truth = set()
    for ui, (c1, c2) in enumerate(pairs):
        h1, h2 = m.units[ui]
        for c, h in ((c1, h1), (c2, h2)):
            for u in adj[c]:
                if u in port_of:
                    truth.add((h, port_of[u]))
    assumptions = [lit if key in truth else -lit for key, lit in m.A.items()]
    r, el, _ = solve_capped(m.cl, 600, assumptions)
    rec["own_assignment_sat"] = r
    rec["own_s"] = round(el, 1)
    # (2) perturbed: swap one port between two connectors in the same block -> expected UNSAT (not
    # guaranteed in general, but a diameter-4 graph is rigid enough here)
    (c0, i0) = min(truth)
    other = next((c, i) for (c, i) in sorted(truth) if c != c0 and PORTS[i][0] == PORTS[i0][0] and (c0, i) not in truth and (c, i0) not in truth)
    pert = (truth - {(c0, i0), other}) | {(c0, other[1]), (other[0], i0)}
    assumptions2 = [lit if key in pert else -lit for key, lit in m.A.items()]
    r2, el2, _ = solve_capped(m.cl, 600, assumptions2)
    rec["perturbed_sat"] = r2
    # (3) free model with symmetry breaking and all redundant clauses: SAT, decoded graph diameter 4
    m2 = Model(len(pairs), [])
    rec["vars"], rec["clauses"] = m2.pool.top, len(m2.cl)
    r3, el3, model = solve_capped(m2.cl, 1800)
    rec["free_sat"], rec["free_s"] = r3, round(el3, 1)
    if r3:
        rec["free_check"] = bfs_check(m2.decode(model))
    # (4) one unit declared as a fixed-type gadget (pure pair): still SAT (exercises gadget code)
    gadget = {"types": [[1, 1, 1], [1, 1, 1]], "edges": [[0, 1]]}
    m3 = Model(len(pairs) - 1, [gadget])
    r4, el4, model4 = solve_capped(m3.cl, 1800)
    rec["gadget_pure_pair_sat"], rec["gadget_s"] = r4, round(el4, 1)
    if r4:
        rec["gadget_check"] = bfs_check(m3.decode(model4))
    print(json.dumps(rec), file=sys.stderr, flush=True)
    json.dump(rec, sys.stdout, indent=1)


def control(args):
    """105 = 15 units + cross triangle {d12, d13, d23}: UNSAT by sat105.py (all 54 skeletons)."""
    tri = {"types": [[1, 1, 0], [1, 0, 1], [0, 1, 1]], "edges": [[0, 1], [0, 2], [1, 2]]}
    cap = float(next((a.split("=")[1] for a in args if a.startswith("--cap=")), 1800))
    # (exact_units, lemma, full, exact_all, tight)
    variants = {"hub": (True, True, False, True, False), "full": (True, True, True, True, False),
                "tight": (True, True, True, True, True), "tight_hub": (True, True, False, True, True),
                "noexact": (True, True, False, False, False), "nounits": (False, True, False, True, False),
                "nolemma": (True, False, False, True, False)}
    names = [a for a in args if a in variants] or ["hub", "full"]
    out = []
    for name in names:
        exact_units, lemma, full, exact_all, tight = variants[name]
        t0 = time.time()
        m = Model(15, [tri], exact_units=exact_units, lemma=lemma, full=full, exact_all=exact_all, tight=tight)
        build = time.time() - t0
        r, el, _ = solve_capped(m.cl, cap)
        rec = {"variant": name, "exact_units": exact_units, "lemma": lemma, "full": full, "exact_all": exact_all, "tight": tight,
               "vars": m.pool.top, "clauses": len(m.cl), "build_s": round(build, 1),
               "result": {True: "SAT", False: "UNSAT", None: "TIMEOUT"}[r], "solve_s": round(el, 1)}
        print(json.dumps(rec), file=sys.stderr, flush=True)
        out.append(rec)
    json.dump(out, sys.stdout, indent=1)


if __name__ == "__main__":
    if sys.argv[1] == "run":
        main_run(sys.argv[2:])
    elif sys.argv[1] == "validate":
        validate(sys.argv[2])
    elif sys.argv[1] == "control":
        control(sys.argv[2:])
