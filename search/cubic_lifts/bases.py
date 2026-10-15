"""Enumerate connected cubic base graphs (multigraphs with loops, parallel edges
and semi-edges) on k vertices up to isomorphism, and write them in the format
read by lift.c:

    k nslots ntree
    u v type          (type 0 = edge u!=v, 1 = loop, 2 = semi-edge; the first
                       ntree slots form a spanning tree and get voltage 1)

A vertex has l loops, s semi-edges and d = 3 - 2l - s ordinary edge ends.
Isomorphism classes are identified by pynauty certificates of a simple graph
encoding: vertex colour = (l, s), each ordinary edge of multiplicity c encoded
by c subdivision vertices.
"""
import itertools
import os
import sys

import pynauty

VTYPES = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1)]  # (loops, semi-edges)


def degree_left(t):
    l, s = t
    return 3 - 2 * l - s


def multigraphs(deg):
    """All symmetric multiplicity matrices (as dict pair->mult) with given degrees."""
    k = len(deg)
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    rem = list(deg)
    cur = {}

    def rec(idx):
        if idx == len(pairs):
            if all(r == 0 for r in rem):
                yield dict(cur)
            return
        i, j = pairs[idx]
        # remaining capacity check: vertex i must be able to place its remaining ends
        # on pairs not yet decided
        maxc = min(rem[i], rem[j])
        for c in range(maxc, -1, -1):
            rem[i] -= c
            rem[j] -= c
            if c:
                cur[(i, j)] = c
            # prune: vertex i has no further pairs after (i, k-1)
            if j == k - 1 and rem[i] != 0:
                pass
            else:
                yield from rec(idx + 1)
            if c:
                del cur[(i, j)]
            rem[i] += c
            rem[j] += c

    yield from rec(0)


def connected(k, mult):
    adj = {i: set() for i in range(k)}
    for (i, j), c in mult.items():
        adj[i].add(j)
        adj[j].add(i)
    seen = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for w in adj[u]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == k


def certificate(k, types, mult):
    nsub = sum(mult.values())
    n = k + nsub
    adj = {i: [] for i in range(n)}
    idx = k
    for (i, j), c in mult.items():
        for _ in range(c):
            adj[i].append(idx)
            adj[j].append(idx)
            idx += 1
    g = pynauty.Graph(n, adjacency_dict=adj)
    colouring = []
    for t in VTYPES:
        cls = {i for i in range(k) if types[i] == t}
        if cls:
            colouring.append(cls)
    if nsub:
        colouring.append(set(range(k, n)))
    g.set_vertex_coloring(colouring)
    return pynauty.certificate(g) + bytes(str(sorted(types)), "ascii")


def spanning_tree(k, mult):
    adj = {i: [] for i in range(k)}
    for (i, j) in mult:
        adj[i].append(j)
        adj[j].append(i)
    seen = {0}
    tree = []
    stack = [0]
    while stack:
        u = stack.pop()
        for w in adj[u]:
            if w not in seen:
                seen.add(w)
                tree.append((min(u, w), max(u, w)))
                stack.append(w)
    return tree


def enumerate_bases(k, allow_semi=True, allow_loops=True):
    out = {}
    tlist = [t for t in VTYPES if (allow_semi or t[1] == 0) and (allow_loops or t[0] == 0)]
    for types in itertools.combinations_with_replacement(tlist, k):
        deg = [degree_left(t) for t in types]
        if sum(deg) % 2:
            continue
        if k > 1 and min(deg) == 0:
            continue
        for mult in multigraphs(deg):
            if k > 1 and not connected(k, mult):
                continue
            c = certificate(k, types, mult)
            if c not in out:
                out[c] = (types, mult)
    return list(out.values())


def write_base(path, k, types, mult):
    tree = spanning_tree(k, mult) if k > 1 else []
    slots = [(u, v, 0) for (u, v) in tree]
    for (u, v), c in sorted(mult.items()):
        extra = c - (1 if (u, v) in tree else 0)
        slots += [(u, v, 0)] * extra
    for i, (l, s) in enumerate(types):
        slots += [(i, i, 1)] * l
        slots += [(i, i, 2)] * s
    with open(path, "w") as f:
        f.write(f"{k} {len(slots)} {len(tree)}\n")
        for u, v, t in slots:
            f.write(f"{u} {v} {t}\n")
    nfree_group = len(slots) - len(tree) - sum(s for _, s in types)
    nsemi = sum(s for _, s in types)
    nloops = sum(l for l, _ in types)
    return nfree_group, nsemi, nloops


if __name__ == "__main__":
    outdir = sys.argv[1]
    ks = [int(x) for x in sys.argv[2:]]
    os.makedirs(outdir, exist_ok=True)
    index = open(os.path.join(outdir, "index.txt"), "a")
    for k in ks:
        bases = enumerate_bases(k)
        nosemi = sum(1 for types, _ in bases if all(s == 0 for _, s in types))
        simple = sum(1 for types, mult in bases if all(t == (0, 0) for t in types) and all(c == 1 for c in mult.values()))
        print(f"k={k}: {len(bases)} bases ({nosemi} without semi-edges, {simple} simple)")
        for i, (types, mult) in enumerate(bases):
            name = f"b{k}_{i:04d}"
            nf, ns, nl = write_base(os.path.join(outdir, name + ".txt"), k, types, mult)
            bundles = sum(1 for c in mult.values() if c > 1)
            index.write(f"{name} k={k} free={nf} semi={ns} loops={nl} bundles={bundles} types={types} mult={sorted(mult.items())}\n")
    index.close()
