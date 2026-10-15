"""Connected typed connector components for three K4,4 blocks with slack, by core + leaves.

Setting as in potential3.py: connector profile x in Z^3, |x| <= 4, connector degree
D[x] = 4 - |x|; edge potential E(x, y) = 2 v(x, y) + node(x) + node(y) <= 0 for all pairs
(checked there), isolated quads have potential self - 15 <= -5. A member on n = 72 + N vertices
has total gadget potential >= B = 96 - 3N and total gadget deficit sum(3 - |x|) = D = 3N - 96
(quads count -1). A component of deficit def therefore needs def - D quads elsewhere when
def > D, each costing >= 5, so its own potential is >= bound(def) = B + 5 max(0, def - D), and
def <= D + floor(D / 5).

Enumeration. Triples have degree 1, so a connected component is (a) two adjacent triples, (b) an
isolated quad, or (c) a connected "core" of doubles / singles / routers (degree 2 / 3 / 4) with
every unfilled core stub carrying a triple leaf. Removing leaves keeps the core connected, and a
triple adjacent to a triple is case (a); so (a)-(c) is a complete case split. Cores are grown
depth-first with a visited set (join two open core vertices or attach a new core vertex),
deduplicated by coloured certificate with (type, core degree) colouring, pruned by deficit <=
maxdef and by the potential bound (E <= 0 everywhere, so the core potential never increases;
singles and routers also pay at least half of max_y E(x, y) < 0 per open stub). Every core state
is completed by all multisets of triple leaves on its open stubs within the potential bound.
Closed components are deduplicated by coloured certificate with type colouring. The old
enumerator (potential3.py) grew triples as internal vertices of partial components and its
frontier exploded at bound -9; cores are bounded by the deficit (at most maxdef vertices), which
is what makes this finite and fast. Its isomorphism key also lacked the colour multiplicities
(see coloured_cert), which dropped 432 components at 106.

Usage: components107.py N [--out=FILE.json] [--flat=BOUND,MAXDEF] [--compare=OLD.json] [--selfcheck]
       (N = number of connectors, 35 for n = 107). shapes_fast.py uses enumerate_components with
       a sink (streaming) and optionally family_bound (sub-family of total potential >= bound).
"""
import collections
import fractions
import hashlib
import itertools
import json
import sys
import time

import pynauty

FR = fractions.Fraction
TYPES = [t for t in itertools.product(range(5), repeat=3) if sum(t) <= 4]
X = {t: sum(t) for t in TYPES}
D = {t: 4 - sum(t) for t in TYPES}
TRIPLES = [t for t in TYPES if X[t] == 3]
CORE = [t for t in TYPES if X[t] <= 2]
QUADS = [t for t in TYPES if X[t] == 4]
NONISO = TRIPLES + CORE


def self_w(x):
    return X[x] ** 2 - sum(v * v for v in x)


def v_w(x, y):
    return X[x] * X[y] - sum(a * b for a, b in zip(x, y))


def node(x):
    return FR(self_w(x) - 3 * X[x] - 3, D[x])


def E(x, y):
    return 2 * v_w(x, y) + node(x) + node(y)


ETAB = {(x, y): E(x, y) for x in NONISO for y in NONISO}
assert max(ETAB.values()) <= 0
BEST = {x: max(ETAB[(x, y)] for y in NONISO) for x in NONISO}
LEAF = {x: sorted(((ETAB[(x, t)], t) for t in TRIPLES), reverse=True) for x in CORE}


def iso_pot(x):
    return FR(self_w(x) - 15)


class Core:
    def __init__(self, types, edges):
        self.types = list(types)
        self.edges = frozenset(edges)
        self.n = len(types)
        self.deg = [0] * self.n
        for a, b in self.edges:
            self.deg[a] += 1
            self.deg[b] += 1
        self.open = [D[self.types[i]] - self.deg[i] for i in range(self.n)]
        self.pot = sum(ETAB[(self.types[a], self.types[b])] for a, b in self.edges)
        self.deficit = sum(3 - X[t] for t in self.types)
        self.slots = tuple(sum(t[j] for t in self.types) for j in range(3))

    def upper(self):
        """Upper bound on the potential of any completion: open stubs of singles and routers cost
        at least half of their best edge each (an edge closes at most two stubs); doubles have
        best 0 (the d-s edge), so only the edges so far count for them."""
        return self.pot + sum(BEST[self.types[i]] * self.open[i] / 2 for i in range(self.n) if X[self.types[i]] <= 1)

    def cert(self):
        return coloured_cert([(t, self.deg[i]) for i, t in enumerate(self.types)], self.edges)


def coloured_cert(colours_of_vertices, edges):
    """Isomorphism key of a vertex-coloured graph. pynauty.certificate is the canonical adjacency
    matrix under a canonical labelling that respects the ordered colour partition, but it does not
    encode the cell sizes: two paths with the same set of colours and different multiplicities can
    share a certificate (found while validating against potential3.py, which had this flaw). The
    key therefore lists every colour with its multiplicity in front of the certificate."""
    colour_of = {}
    for i, c in enumerate(colours_of_vertices):
        colour_of.setdefault(c, set()).add(i)
    keys = sorted(colour_of)
    adj = {i: [] for i in range(len(colours_of_vertices))}
    for a, b in edges:
        adj[a].append(b)
    g = pynauty.Graph(len(colours_of_vertices), directed=False, adjacency_dict=adj, vertex_coloring=[colour_of[k] for k in keys])
    return (tuple((k, len(colour_of[k])) for k in keys), pynauty.certificate(g))


def closed_cert(types, edges):
    return coloured_cert(list(types), edges)


def leaf_completions(core, bound):
    """All ways to put triple leaves on the open stubs of core with total potential >= bound."""
    stubs = [(i, core.open[i]) for i in range(core.n) if core.open[i] > 0]
    budget = core.pot - bound  # total leaf cost allowed (>= 0 required)
    if budget < 0:
        return
    # per-vertex multisets of leaf types with their cost, cheapest first
    per_vertex = []
    for i, o in stubs:
        opts = []
        for combo in itertools.combinations_with_replacement(LEAF[core.types[i]], o):
            cost = -sum(e for e, _ in combo)
            if cost <= budget:
                opts.append((cost, tuple(t for _, t in combo)))
        opts.sort()
        per_vertex.append(opts)
    # minimum cost of the remaining vertices, for pruning
    mins = [opts[0][0] if opts else None for opts in per_vertex]
    if any(m is None for m in mins):
        return
    suffix_min = [0] * (len(stubs) + 1)
    for k in range(len(stubs) - 1, -1, -1):
        suffix_min[k] = suffix_min[k + 1] + mins[k]

    def rec(k, spent, chosen):
        if k == len(stubs):
            yield list(chosen)
            return
        for cost, leaves in per_vertex[k]:
            if spent + cost + suffix_min[k + 1] > budget:
                break
            chosen.append(leaves)
            yield from rec(k + 1, spent + cost, chosen)
            chosen.pop()

    for choice in rec(0, 0, []):
        types = list(core.types)
        edges = set(core.edges)
        for (i, _), leaves in zip(stubs, choice):
            for t in leaves:
                edges.add((i, len(types)))
                types.append(t)
        yield types, edges


def enumerate_components(N, log=sys.stderr, flat=None, sink=None, family_bound=None):
    """flat = (bound, maxdef): constant bound for every deficit (to reproduce potential3.py runs).
    sink(types, edges, potential): streaming mode; components are passed on as they are found
    (deduplicated by certificate, stored as 16-byte digests of the key) and only the summary is
    returned. family_bound: restrict to shapes of total potential >= family_bound (> B): the
    component bound becomes family_bound + 5 max(0, d - D), a sub-family of the full problem."""
    Dtot = 3 * N - 96
    B = FR(96 - 3 * N)
    if family_bound is not None:
        B = FR(family_bound)
    maxdef = Dtot + int(-B) // 5
    if flat:
        B, maxdef = FR(flat[0]), flat[1]

    def bound(deficit):
        if flat:
            return B
        return B + 5 * max(0, deficit - Dtot)

    found = {} if sink is None else set()
    stats = collections.Counter()
    by_pot = collections.Counter()
    by_def = collections.Counter()

    def emit(types, edges, pot, kind):
        k = closed_cert(types, edges)
        if sink is not None:
            k = hashlib.blake2b(repr(k).encode(), digest_size=16).digest()
        if k in found:
            stats["dup_" + kind] += 1
            return
        stats[kind] += 1
        if sink is None:
            found[k] = {"types": [list(t) for t in types], "edges": sorted(edges), "potential": str(pot), "n": len(types),
                        "slots": [sum(t[j] for t in types) for j in range(3)], "kind": kind}
        else:
            found.add(k)
            by_pot[str(pot)] += 1
            by_def[sum(3 - X[t] for t in types)] += 1
            sink(list(types), set(edges), pot)

    # (a) triple pairs, (b) quads
    for a, b in itertools.combinations_with_replacement(TRIPLES, 2):
        if ETAB[(a, b)] >= bound(0):
            emit([a, b], {(0, 1)}, ETAB[(a, b)], "pair")
    for q in QUADS:
        if iso_pot(q) >= bound(-1):
            emit([q], set(), iso_pot(q), "quad")
    # (c) cores with leaves. Depth-first over core states with a visited set (digests of the
    # coloured certificate); every state is expanded exactly once, as in the breadth-first version,
    # but only the stack of pending states is held in memory (the breadth-first frontier at 108
    # exceeded 2 GB at level 9).
    def state_key(c):
        return hashlib.blake2b(repr(c.cert()).encode(), digest_size=16).digest()

    seen = set()
    stack = []
    for t in CORE:
        c = Core([t], [])
        seen.add(state_key(c))
        stack.append(c)
    t0 = time.time()
    expanded = 0
    next_report = 100000
    while stack:
        c = stack.pop()
        expanded += 1
        # complete this core with leaves
        for types, edges in leaf_completions(c, bound(c.deficit)):
            pot = sum(ETAB[(types[a], types[b])] for a, b in edges)
            assert pot >= bound(c.deficit)
            emit(types, edges, pot, "core")
        opens = [i for i in range(c.n) if c.open[i] > 0]
        cands = []
        for a, b in itertools.combinations(opens, 2):
            if (a, b) not in c.edges:
                cands.append((c.types, c.edges | {(a, b)}))
        if c.deficit < maxdef:
            for a in opens:
                for t in CORE:
                    if c.deficit + 3 - X[t] > maxdef:
                        continue
                    if any(c.slots[j] + t[j] > 32 for j in range(3)):
                        continue
                    cands.append((c.types + [t], c.edges | {(a, c.n)}))
        for types, edges in cands:
            c2 = Core(types, edges)
            if c2.upper() < bound(c2.deficit):
                stats["pruned_potential"] += 1
                continue
            k = state_key(c2)
            if k in seen:
                stats["dup_core_state"] += 1
                continue
            seen.add(k)
            stack.append(c2)
        if expanded >= next_report:
            next_report += 100000
            print(f"expanded {expanded} core states, stack {len(stack)}, {len(seen)} seen, {len(found)} components, "
                  f"{time.time() - t0:.0f} s", file=log, flush=True)
    print(f"done: {expanded} core states, {len(found)} components, {time.time() - t0:.0f} s", file=log, flush=True)
    if sink is not None:
        return {"N": N, "n": 72 + N, "D": Dtot, "B": str(B), "maxdef": maxdef, "components": len(found),
                "core_states": len(seen), "stats": dict(stats),
                "by_potential": {k: by_pot[k] for k in sorted(by_pot, key=FR, reverse=True)},
                "by_deficit": {str(k): by_def[k] for k in sorted(by_def)}, "wall_s": round(time.time() - t0, 1)}
    comps = list(found.values())
    comps.sort(key=lambda c: (-FR(c["potential"]), c["n"], c["types"]))
    by_pot = collections.Counter(c["potential"] for c in comps)
    by_def = collections.Counter(sum(3 - sum(t) for t in c["types"]) for c in comps)
    summary = {"N": N, "n": 72 + N, "D": Dtot, "B": str(B), "maxdef": maxdef, "components": len(comps),
               "core_states": len(seen), "stats": dict(stats),
               "by_potential": {k: by_pot[k] for k in sorted(by_pot, key=FR, reverse=True)},
               "by_deficit": {str(k): by_def[k] for k in sorted(by_def)},
               "by_kind": dict(collections.Counter(c["kind"] for c in comps))}
    return comps, summary


def selfcheck(comps, max_group=400):
    """Independent check of the certificate dedup with networkx: within every type multiset the
    listed components must be pairwise non-isomorphic as typed graphs (groups above max_group
    members are skipped and reported)."""
    import networkx as nx

    def G(c):
        g = nx.Graph()
        for i, t in enumerate(c["types"]):
            g.add_node(i, t=tuple(t))
        g.add_edges_from(tuple(e) for e in c["edges"])
        return g

    nm = nx.algorithms.isomorphism.categorical_node_match("t", None)
    groups = collections.defaultdict(list)
    for c in comps:
        groups[tuple(sorted(tuple(t) for t in c["types"]))].append(c)
    dup, tested, skipped = 0, 0, 0
    for key, grp in groups.items():
        if len(grp) > max_group:
            skipped += 1
            continue
        gs = [G(c) for c in grp]
        for i in range(len(gs)):
            for j in range(i + 1, len(gs)):
                tested += 1
                if nx.is_isomorphic(gs[i], gs[j], node_match=nm):
                    dup += 1
    return {"groups": len(groups), "pairs_tested": tested, "isomorphic_pairs": dup, "groups_skipped": skipped}


def compare_with(path, comps):
    """Certificate-set comparison with a components file of potential3.py (validation)."""
    old = json.load(open(path))
    oldc = {closed_cert([tuple(t) for t in c["types"]], {tuple(e) for e in c["edges"]}): c for c in old}
    newc = {closed_cert([tuple(t) for t in c["types"]], {tuple(e) for e in c["edges"]}): c for c in comps}
    only_old = [oldc[k] for k in oldc.keys() - newc.keys()]
    only_new = [newc[k] for k in newc.keys() - oldc.keys()]
    return {"old": len(oldc), "new": len(newc), "common": len(oldc.keys() & newc.keys()),
            "only_old": only_old[:5], "only_new": only_new[:5], "n_only_old": len(only_old), "n_only_new": len(only_new)}


if __name__ == "__main__":
    N = int(sys.argv[1])
    out = next((a.split("=")[1] for a in sys.argv if a.startswith("--out=")), None)
    flat = next((a.split("=")[1] for a in sys.argv if a.startswith("--flat=")), None)
    flat = tuple(int(x) for x in flat.split(",")) if flat else None
    cmp = next((a.split("=")[1] for a in sys.argv if a.startswith("--compare=")), None)
    t0 = time.time()
    comps, summary = enumerate_components(N, flat=flat)
    summary["wall_s"] = round(time.time() - t0, 1)
    if cmp:
        summary["compare"] = compare_with(cmp, comps)
    if "--selfcheck" in sys.argv:
        summary["selfcheck"] = selfcheck(comps)
    print(json.dumps(summary, indent=1), file=sys.stderr, flush=True)
    if out:
        with open(out, "w") as f:
            json.dump({"summary": summary, "components": comps}, f)
    else:
        json.dump({"summary": summary, "components": comps}, sys.stdout)
