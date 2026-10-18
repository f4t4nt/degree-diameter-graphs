#!/usr/bin/env python3
"""Controllers as regular lifts of a k-vertex Schreier base over a group G (multiplication table).
Symbols: 0 = a (dart c -> alpha0[c], voltage va[c]), 1 = a^-1, 2 = b (c -> beta0[c], vb[c]), 3 = b^-1,
         4 = s (c -> sigma0[c], vs[c]; sigma0 an involution of the base: c <-> sigma0[c] is a link with voltage
         vs[c] on the dart from the smaller index and its inverse on the reverse dart; a fixed point c of sigma0
         is a semi-edge whose voltage vs[c] must be an involution of G).
Lift convention: dart c -> d with voltage v joins (c, g) to (d, g v); the deck group acts on the left, so the
k vertices (c, e) are root representatives.  Vertex (c, g) has index c * |G| + g (identity e = 0).
"""
import json

def cyclic_table(n):
    return [[(a + b) % n for b in range(n)] for a in range(n)]

def read_table(fn):
    with open(fn) as f:
        toks = f.read().split()
    n = int(toks[0])
    vals = list(map(int, toks[1:1 + n * n]))
    return [vals[i * n:(i + 1) * n] for i in range(n)]

def inverses(T):
    n = len(T)
    return [next(h for h in range(n) if T[g][h] == 0) for g in range(n)]

class Base:
    def __init__(self, k, alpha0, beta0, sigma0):
        self.k = k; self.alpha0 = list(alpha0); self.beta0 = list(beta0)
        self.sigma0 = list(sigma0) if sigma0 is not None else None
        if self.sigma0 is not None:
            assert all(self.sigma0[self.sigma0[c]] == c for c in range(k)), 'sigma0 must be an involution'

def lift_controller(base, T, va, vb, vs, symbols=(1, 0, 3, 2, 4)):
    """returns N, out[v] = [(target, symbol)...].  symbols: the alphabet's inv list, used only to decide whether
    the s symbol (index 4) is present."""
    nG = len(T); inv = inverses(T); k = base.k; N = k * nG
    a0, b0, s0 = base.alpha0, base.beta0, base.sigma0
    use_s = len(symbols) >= 5
    alpha = [0] * N; beta = [0] * N; sig = [0] * N
    for c in range(k):
        for g in range(nG):
            v = c * nG + g
            alpha[v] = a0[c] * nG + T[g][va[c]]
            beta[v] = b0[c] * nG + T[g][vb[c]]
            if use_s:
                d = s0[c]
                if d == c:
                    assert T[vs[c]][vs[c]] == 0 and vs[c] != 0, 'semi-edge voltage must be a non-trivial involution'
                    sig[v] = c * nG + T[g][vs[c]]
                elif c < d:
                    sig[v] = d * nG + T[g][vs[c]]
                else:
                    sig[v] = d * nG + T[g][inv[vs[d]]]
    ainv = [0] * N; binv = [0] * N
    for v in range(N): ainv[alpha[v]] = v; binv[beta[v]] = v
    if use_s:
        assert all(sig[sig[v]] == v for v in range(N))
    out = []
    for v in range(N):
        o = [(alpha[v], 0), (ainv[v], 1), (beta[v], 2), (binv[v], 3)]
        if use_s: o.append((sig[v], 4))
        out.append(o)
    return N, out

def controller_stats(N, out):
    """simple?, degree, diameter, girth (<= 5 detection), distance distribution from vertex 0."""
    adj = [[t for t, _ in o] for o in out]
    simple = all(len(set(a)) == len(a) and v not in a for v, a in enumerate(adj))
    from collections import deque
    def bfs(s):
        dist = [-1] * N; dist[s] = 0; q = deque([s])
        while q:
            x = q.popleft()
            for y in adj[x]:
                if dist[y] < 0: dist[y] = dist[x] + 1; q.append(y)
        return dist
    d0 = bfs(0)
    diam = max(max(bfs(s)) for s in range(N))
    hist = [d0.count(i) for i in range(max(d0) + 1)]
    return simple, len(adj[0]), diam, hist

def save_controller_json(fn, base, T_name, va, vb, vs, N, out):
    json.dump({'k': base.k, 'alpha0': base.alpha0, 'beta0': base.beta0, 'sigma0': base.sigma0, 'group': T_name,
               'va': va, 'vb': vb, 'vs': vs, 'N': N, 'edges': [[c, d, s] for c in range(N) for d, s in out[c]]}, open(fn, 'w'))

if __name__ == '__main__':
    # write a controller as an exlift.c controller file: N K / inv / dim (all line symbols) / nbr rows / the k roots (c, e)
    import sys
    if len(sys.argv) < 10:
        sys.exit('usage: ctrl.py group.tbl k alpha0 beta0 sigma0|- va vb vs|- out.ctrl')
    T = read_table(sys.argv[1]); k = int(sys.argv[2])
    a0, b0 = (list(map(int, s.split(','))) for s in sys.argv[3:5])
    s0 = None if sys.argv[5] == '-' else list(map(int, sys.argv[5].split(',')))
    va, vb = (list(map(int, s.split(','))) for s in sys.argv[6:8])
    vs = None if sys.argv[8] == '-' else list(map(int, sys.argv[8].split(',')))
    N, out = lift_controller(Base(k, a0, b0, s0), T, va, vb, vs, symbols=(1, 0, 3, 2, 4) if s0 is not None else (1, 0, 3, 2))
    K = len(out[0])
    with open(sys.argv[9], 'w') as f:
        f.write('%d %d\n' % (N, K))
        f.write(' '.join(map(str, [1, 0, 3, 2, 4][:K])) + '\n')
        f.write(' '.join(['1'] * K) + '\n')
        for o in out:
            f.write(' '.join(str(next(d for d, s in o if s == t)) for t in range(K)) + '\n')
        f.write('%d%s\n' % (k, ''.join(' %d' % (c * len(T)) for c in range(k))))
    print('controller (simple, degree, diameter, distances from 0):', controller_stats(N, out))
