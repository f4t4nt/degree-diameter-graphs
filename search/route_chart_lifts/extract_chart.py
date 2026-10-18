#!/usr/bin/env python3
"""Recover the controller + route chart of a route-chart lift from its edge list.

Vertex label convention of the published files: v = F*c + x, c = controller vertex,
x in [0, F) the fibre coordinate written in base p (digit i = coordinate i), F = p^m.
For each ordered pair of adjacent fibres (c, d) the neighbour set of (c, x) in fibre d
must be an affine subspace  A x + T  with A linear (mod T) and T a fixed r-dim subspace.
Output (JSON): controller edges with chart ids, charts (A as m x m matrix over F_p, T basis).
Usage: extract_chart.py edges p m out.json
"""
import sys, json, itertools, collections

def digits(x, p, m):
    d = []
    for _ in range(m):
        d.append(x % p); x //= p
    return d

def undigits(d, p):
    x = 0
    for i in reversed(range(len(d))):
        x = x * p + d[i]
    return x

def vadd(a, b, p):
    return [(u + v) % p for u, v in zip(a, b)]

def vsub(a, b, p):
    return [(u - v) % p for u, v in zip(a, b)]

def vscale(l, a, p):
    return [(l * u) % p for u in a]

def span(vecs, p, m):
    S = {tuple([0] * m)}
    for v in vecs:
        new = set()
        for s in S:
            for l in range(p):
                new.add(tuple(vadd(list(s), vscale(l, v, p), p)))
        S = new
    return S

def rank(vecs, p, m):
    # gaussian elimination over F_p
    rows = [list(v) for v in vecs]
    r = 0
    for col in range(m):
        piv = None
        for i in range(r, len(rows)):
            if rows[i][col] % p:
                piv = i; break
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][col], -1, p)
        rows[r] = [(inv * u) % p for u in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][col] % p:
                f = rows[i][col]
                rows[i] = [(u - f * w) % p for u, w in zip(rows[i], rows[r])]
        r += 1
    return r

def basis_of_subspace(S, p, m):
    vecs = []
    for v in sorted(S):
        if any(v) and rank(vecs + [list(v)], p, m) > len(vecs):
            vecs.append(list(v))
    return vecs

def main():
    fn, p, m, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    F = p ** m
    adj = collections.defaultdict(set)
    for line in open(fn):
        if not line.strip() or line.startswith('#'): continue
        a, b = map(int, line.split()); adj[a].add(b); adj[b].add(a)
    n = len(adj); N = n // F
    assert n == N * F, (n, F)
    # controller adjacency: fibre pairs
    cadj = collections.defaultdict(lambda: collections.defaultdict(set))
    for v in range(n):
        for w in adj[v]:
            cadj[v // F][w // F].add((v % F, w % F))
    print('n =', n, 'N =', N, 'F =', F)
    degs = collections.Counter(len(cadj[c]) for c in range(N))
    print('controller degree histogram:', dict(degs))
    charts = []   # list of (A rows, T basis) canonical
    chart_index = {}
    edges = []    # (c, d, chart id)
    zero = tuple([0] * m)
    for c in range(N):
        for d in sorted(cadj[c]):
            pairs = cadj[c][d]
            Y = collections.defaultdict(set)
            for x, y in pairs: Y[x].add(y)
            sizes = {len(s) for s in Y.values()}
            assert len(Y) == F and len(sizes) == 1, (c, d, sizes, len(Y))
            # translation subspace from x = 0
            Y0 = sorted(Y[0]); y0 = Y0[0]
            T = {tuple(vsub(digits(y, p, m), digits(y0, p, m), p)) for y in Y0}
            Tb = basis_of_subspace(T, p, m)
            assert span(Tb, p, m) == T, 'Y_0 not an affine subspace'
            r = len(Tb)
            # choose A e_i in Y[e_i] (any representative), check A x + T == Y[x] for all x
            Acols = []
            for i in range(m):
                e = [0] * m; e[i] = 1
                yi = sorted(Y[undigits(e, p)])[0]
                Acols.append(digits(yi, p, m))
            ok = True
            for x in range(F):
                xd = digits(x, p, m)
                Ax = [0] * m
                for i in range(m):
                    if xd[i]: Ax = vadd(Ax, vscale(xd[i], Acols[i], p), p)
                img = {undigits(vadd(Ax, list(t), p), p) for t in T}
                if img != Y[x]: ok = False; break
            assert ok, ('not affine-linear', c, d)
            # canonical form: reduce columns of A modulo T so equal charts compare equal:
            # for each column pick the lexicographically smallest element of col + T
            Acan = []
            for col in Acols:
                Acan.append(min(tuple(vadd(col, list(t), p)) for t in T))
            key = (tuple(Acan), tuple(sorted(T)))
            if key not in chart_index:
                chart_index[key] = len(charts)
                charts.append({'A_cols_canonical': [list(c_) for c_ in Acan], 'A_cols': Acols, 'T': sorted(list(t) for t in T), 'T_basis': Tb, 'r': r})
            edges.append((c, d, chart_index[key]))
    print('directed controller edges:', len(edges), ' distinct charts:', len(charts))
    cnt = collections.Counter(e[2] for e in edges)
    print('chart usage:', sorted(cnt.items()))
    # inverse pairing: chart of (d,c) vs (c,d)
    emap = {(c, d): k for c, d, k in edges}
    inv = collections.Counter()
    for (c, d), k in emap.items():
        inv[(k, emap[(d, c)])] += 1
    print('chart inverse pairs (k, k^-1): count', sorted(inv.items()))
    json.dump({'p': p, 'm': m, 'N': N, 'edges': edges, 'charts': charts}, open(out, 'w'))
    print('wrote', out)

main()
