"""Independent verifier and three-block decomposition for degree-4 diameter-4 candidates.

Usage:
  python verify.py corpus            # the 36 diameter-4 graphs on 104 vertices in ../classes/
  python verify.py edges FILE...     # edge-list files ("u v" per line, '#' lines ignored)

For every graph: degree/simplicity/connectivity/all-pairs diameter by BFS, then a search for three
vertex-disjoint subdivided K4,4 blocks covering 72 vertices. For a found decomposition it reports the
connector port profiles, the connector graph, and checks the hub-to-port reach identity

    sum_{hub a} sum_{block j != block(a)} S_j(a)  ==  sum_c self(c) + sum_{edges cc'} w(c, c')

where S_j(a) = sum over the 8 slot connectors c of a's ports of [x_j(c) + sum_{c' ~ c} x_j(c')],
self(c) = 2 (x(c)^2 - sum_j x_j(c)^2) and w(c,c') = 2 sum_j [(x(c)-x_j(c)) x_j(c') + (x(c')-x_j(c')) x_j(c)].
"""
import collections
import itertools
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
CLASSES = HERE.parent / "classes"


def load_corpus():
    """The 36 isomorphism classes at 104 (classes/classNN.edges.txt) in the shape the campaign's
    corpus file had: {"representatives": [{"edges_sorted": [[u, v], ...]}, ...]}, index = class."""
    reps = []
    for path in sorted(CLASSES.glob("class*.edges.txt")):
        reps.append({"name": path.stem.split(".")[0], "edges_sorted": [list(e) for e in load_edges_file(path)]})
    return {"representatives": reps}


def load_edges_file(path):
    edges = []
    for line in pathlib.Path(path).read_text().split("\n"):
        parts = line.split()
        if len(parts) >= 2 and parts[0].lstrip("-").isdigit() and parts[1].lstrip("-").isdigit():
            edges.append((int(parts[0]), int(parts[1])))
    return edges


def build(edges):
    adj = collections.defaultdict(set)
    seen = set()
    problems = []
    for u, v in edges:
        if u == v:
            problems.append(f"self-loop {u}")
        key = (min(u, v), max(u, v))
        if key in seen:
            problems.append(f"duplicate edge {key}")
        seen.add(key)
        adj[u].add(v)
        adj[v].add(u)
    verts = sorted(adj)
    return verts, adj, problems


def bfs(adj, src):
    dist = {src: 0}
    frontier = [src]
    while frontier:
        nxt = []
        for u in frontier:
            for v in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + 1
                    nxt.append(v)
        frontier = nxt
    return dist


def basic_checks(verts, adj):
    n = len(verts)
    degs = collections.Counter(len(adj[v]) for v in verts)
    bad_pairs = 0
    diameter = 0
    connected = True
    for v in verts:
        d = bfs(adj, v)
        if len(d) != n:
            connected = False
        far = max(d.values())
        diameter = max(diameter, far)
        bad_pairs += sum(1 for u, x in d.items() if x > 4 and u > v)
    return {"n": n, "degrees": dict(degs), "connected": connected, "diameter": diameter,
            "unordered_pairs_beyond_4": bad_pairs}


def k44_blocks(verts, adj):
    """All subdivided-K4,4 blocks: frozenset of 24 vertices -> (L hubs, R hubs, ports dict)."""
    found = {}
    for l0 in verts:
        if len(adj[l0]) != 4:
            continue
        ports0 = sorted(adj[l0])
        choices = [sorted(adj[p] - {l0}) for p in ports0]
        for rs in itertools.product(*choices):
            if len(set(rs)) != 4:
                continue
            R = set(rs)
            # every R hub must see 4 ports whose other hub endpoints form one common L set of size 4
            L = None
            ok = True
            port_of = {}
            for r in R:
                if len(adj[r]) != 4:
                    ok = False
                    break
                lset = set()
                for q in adj[r]:
                    others = adj[q] - {r}
                    # q must be adjacent to exactly one hub of L; L unknown yet, collect candidates
                    lset_q = others
                    port_of[(r, q)] = lset_q
                # defer; collect union
            if not ok:
                continue
            # candidate L: l0 plus hubs that are, for each r, endpoints of r's ports
            cand = set.intersection(*[set().union(*[adj[q] - {r} for q in adj[r]]) for r in R])
            cand.add(l0)
            if len(cand) < 4:
                continue
            for L in itertools.combinations(sorted(cand), 4):
                if l0 not in L:
                    continue
                Lset = set(L)
                if Lset & R:
                    continue
                ports = {}
                good = True
                for l in L:
                    if len(adj[l]) != 4:
                        good = False
                        break
                    for q in adj[l]:
                        rr = adj[q] & R
                        if len(rr) != 1 or len(adj[q] & Lset) != 1:
                            good = False
                            break
                        ports[(l, next(iter(rr)))] = q
                    if not good:
                        break
                if not good or len(ports) != 16:
                    continue
                if any(len(adj[r] - set(ports.values())) for r in R):
                    continue
                block = frozenset(Lset | R | set(ports.values()))
                if len(block) == 24 and block not in found:
                    found[block] = (tuple(sorted(Lset)), tuple(sorted(R)), {k: v for k, v in ports.items()})
    return found


def decompositions(blocks, n):
    keys = list(blocks)
    out = []
    for trip in itertools.combinations(keys, 3):
        if len(trip[0] | trip[1] | trip[2]) == 72:
            out.append(trip)
    return out


def analyze_decomposition(verts, adj, trip, blocks):
    block_of = {}
    hubs, ports = {}, {}
    for i, b in enumerate(trip):
        L, R, pmap = blocks[b]
        for v in L + R:
            block_of[v] = i
            hubs[v] = i
        for q in pmap.values():
            block_of[q] = i
            ports[q] = i
    connectors = [v for v in verts if v not in block_of]
    cset = set(connectors)
    x = {c: [0, 0, 0] for c in connectors}
    cadj = {c: set() for c in connectors}
    bad = []
    for c in connectors:
        for u in adj[c]:
            if u in ports:
                x[c][ports[u]] += 1
            elif u in cset:
                cadj[c].add(u)
            else:
                bad.append((c, u))
    profiles = collections.Counter()
    for c in connectors:
        profiles[(sum(x[c]), tuple(sorted(x[c], reverse=True)))] += 1
    # reach identity
    lhs = 0
    reach_ok = True
    for a, i in hubs.items():
        slots = []
        for q in adj[a]:
            slots.extend(sorted(adj[q] & cset))
        for j in range(3):
            if j == i:
                continue
            multiset = collections.Counter()
            for c in slots:
                for q in adj[c]:
                    if q in ports and ports[q] == j:
                        multiset[q] += 1
                for c2 in cadj[c]:
                    for q in adj[c2]:
                        if q in ports and ports[q] == j:
                            multiset[q] += 1
            lhs += sum(multiset.values())
            if len(multiset) != 16:
                reach_ok = False
    rhs = 0
    for c in connectors:
        xc = sum(x[c])
        rhs += 2 * (xc * xc - sum(t * t for t in x[c]))
    for c in connectors:
        for c2 in cadj[c]:
            if c < c2:
                xc, xc2 = sum(x[c]), sum(x[c2])
                rhs += 2 * sum((xc - x[c][j]) * x[c2][j] + (xc2 - x[c2][j]) * x[c][j] for j in range(3))
    return {"connectors": len(connectors), "non_port_non_connector_edges": len(bad),
            "profiles": {f"x={k[0]} {k[1]}": v for k, v in sorted(profiles.items())},
            "connector_degree_hist": dict(collections.Counter(len(cadj[c]) for c in connectors)),
            "reach_identity_lhs": lhs, "reach_identity_rhs": rhs,
            "every_hub_reaches_all_16_ports_of_each_other_block": reach_ok}


def run_graph(name, edges):
    verts, adj, problems = build(edges)
    rec = {"name": name, "edges": len(edges), "problems": problems}
    rec.update(basic_checks(verts, adj))
    blocks = k44_blocks(verts, adj)
    decs = decompositions(blocks, len(verts))
    rec["k44_blocks_found"] = len(blocks)
    rec["three_block_decompositions"] = len(decs)
    if decs:
        rec["decomposition"] = analyze_decomposition(verts, adj, decs[0], blocks)
        rec["decomposition_block_sets"] = [sorted(b) for b in decs[0]]
    return rec


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "corpus"
    results = []
    if mode == "corpus":
        corpus = load_corpus()
        for i, rep in enumerate(corpus["representatives"]):
            results.append(run_graph(f"class{i}", [tuple(e) for e in rep["edges_sorted"]]))
    else:
        for path in sys.argv[2:]:
            results.append(run_graph(path, load_edges_file(path)))
    print(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
