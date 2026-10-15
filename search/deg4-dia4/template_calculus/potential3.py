"""Potential function and connector-component enumeration for the three-K4,4-block family with slack.

Connector type = port profile x in Z^3 (ports per block), |x| <= 4, connector degree d = 4 - |x|
(4-regular, all 96 slots filled). F(c) = x(c) + sum of x over connector neighbours,
T(c) = |x||F| - x.F = self(x) + sum_{c' ~ c} v(x, x'), self(x) = |x|^2 - x.x, v(x, y) = |x||y| - x.y.
Hub-to-port routes (TEMPLATE_CALCULUS.md section 4): sum over the 48 (hub, foreign block) pairs of
S_j(a) equals 2 sum_c T(c), and every S_j(a) >= 16, so sum_c T(c) >= 384.

phi(c) = T(c) - 3|x| - 3. With N connectors and 96 slots, sum_c phi(c) = sum_c T(c) - 288 - 3N, so a
member on 72 + N vertices has sum_c phi(c) >= 96 - 3N (-6 at n = 106, -9 at 107, -12 at 108).
Distributing node potentials over edges, node(c) = (self - 3|x| - 3)/d, E(c, c') = 2 v + node(c) +
node(c'), gives sum_c phi = sum_edges E + sum_isolated (self - 15). This script tabulates E (all
E <= 0 is checked, not assumed) and enumerates all connected typed connector graphs ("components",
every vertex saturated: degree 4 - |x|) with potential >= bound, up to isomorphism (pynauty
certificates with type colouring).

Usage: potential3.py table | components BOUND MAXDEF
"""
import fractions
import itertools
import json
import sys

import pynauty

TYPES = [t for t in itertools.product(range(5), repeat=3) if sum(t) <= 4]
X = {t: sum(t) for t in TYPES}
D = {t: 4 - sum(t) for t in TYPES}
FR = fractions.Fraction


def self_w(x):
    return X[x] ** 2 - sum(v * v for v in x)


def v_w(x, y):
    return X[x] * X[y] - sum(a * b for a, b in zip(x, y))


def node(x):
    return FR(self_w(x) - 3 * X[x] - 3, D[x])


def E(x, y):
    return 2 * v_w(x, y) + node(x) + node(y)


def iso_pot(x):
    return self_w(x) - 15


NONISO = [t for t in TYPES if D[t] > 0]
ETAB = {(x, y): E(x, y) for x in NONISO for y in NONISO}


def table():
    out = {"types": len(TYPES), "nonisolated": len(NONISO)}
    worst = max(ETAB.values())
    out["max_E"] = str(worst)
    out["all_E_nonpositive"] = worst <= 0
    zero = sorted({tuple(sorted((x, y))) for (x, y), e in ETAB.items() if e == 0})
    out["E_zero_pairs"] = [f"{x}-{y}" for x, y in zero]
    out["isolated_potentials"] = {str(t): iso_pot(t) for t in TYPES if D[t] == 0}
    small = sorted({(str(e), f"{x}-{y}") for (x, y), e in ETAB.items() if x <= y and -6 <= e < 0})
    out["E_in_[-6,0)"] = small
    print(json.dumps(out, indent=1))


class Comp:
    """A partial typed connector graph: types[i], adjacency sets, open degrees."""

    def __init__(self, types, edges):
        self.types = list(types)
        self.edges = set(edges)
        self.n = len(self.types)
        self.deg = [0] * self.n
        for a, b in self.edges:
            self.deg[a] += 1
            self.deg[b] += 1
        self.open = [D[self.types[i]] - self.deg[i] for i in range(self.n)]
        assert all(o >= 0 for o in self.open)

    def potential(self):
        if self.n == 1 and D[self.types[0]] == 0:
            return FR(iso_pot(self.types[0]))
        return sum(ETAB[(self.types[a], self.types[b])] for a, b in self.edges)

    def slots(self):
        return tuple(sum(t[j] for t in self.types) for j in range(3))

    def closed(self):
        return all(o == 0 for o in self.open)

    def cert(self):
        colour_of = {}
        for i, t in enumerate(self.types):
            colour_of.setdefault((t, self.open[i]), set()).add(i)
        colours = [colour_of[k] for k in sorted(colour_of)]
        adj = {i: [] for i in range(self.n)}
        for a, b in self.edges:
            adj[a].append(b)
        g = pynauty.Graph(self.n, directed=False, adjacency_dict=adj, vertex_coloring=colours)
        return (tuple(sorted(colour_of)), pynauty.certificate(g))


def components(bound, maxdef):
    """All closed connected components with potential >= bound and slot deficit sum(3 - |x|) <= maxdef
    (with N connectors on 96 slots the total deficit is 3N - 96: 6 at n = 106, 9 at 107, 12 at 108;
    quads have deficit -1 and potential -5, so a component's deficit is at most the total plus the
    number of quads the potential allows). Breadth-first growth of partial components (deduplicated
    by certificate): either join two open vertices by a new edge or attach a new vertex to an open
    vertex. Non-isolated types have |x| <= 3, so the deficit of a partial component never decreases."""
    bound = FR(bound)
    seen = set()
    closed = {}
    frontier = []
    for t in TYPES:
        c = Comp([t], [])
        k = c.cert()
        seen.add(k)
        if c.closed():
            if c.potential() >= bound:
                closed[k] = c
        else:
            frontier.append(c)
    while frontier:
        nxt = []
        for c in frontier:
            opens = [i for i in range(c.n) if c.open[i] > 0]
            cands = []
            for a, b in itertools.combinations(opens, 2):
                if (a, b) not in c.edges:
                    cands.append((c.types, c.edges | {(a, b)}))
            cdef = sum(3 - X[t] for t in c.types)
            for a in opens:
                for t in NONISO:
                    if cdef + 3 - X[t] > maxdef:
                        continue
                    sl = c.slots()
                    if any(sl[j] + t[j] > 32 for j in range(3)):
                        continue
                    cands.append((c.types + [t], c.edges | {(a, c.n)}))
            for types, edges in cands:
                c2 = Comp(types, edges)
                if c2.potential() < bound:
                    continue
                k = c2.cert()
                if k in seen:
                    continue
                seen.add(k)
                if c2.closed():
                    closed[k] = c2
                else:
                    nxt.append(c2)
        frontier = nxt
        print(f"frontier {len(frontier)} partial, {len(closed)} closed, {len(seen)} seen", file=sys.stderr, flush=True)
    return list(closed.values())


def comp_record(c):
    return {"types": [list(t) for t in c.types], "edges": sorted(c.edges), "potential": str(c.potential()),
            "n": c.n, "slots": c.slots()}


if __name__ == "__main__":
    if sys.argv[1] == "table":
        table()
    elif sys.argv[1] == "components":
        bound = FR(sys.argv[2])
        maxdef = int(sys.argv[3])
        comps = components(bound, maxdef)
        comps.sort(key=lambda c: (-c.potential(), c.n, c.types))
        print(json.dumps([comp_record(c) for c in comps], indent=0))
