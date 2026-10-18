#!/usr/bin/env python3
"""General route-chart alphabets: inverse pairs (t, t^-1) AND self-inverse (involution) symbols.

Symbol t carries a chart (A_t, d_t): A_t in GL_m(F_p) (columns), translation line T_t = F_p d_t.
Lift rule: (v, x) ~ (t(v), A_t x + lambda d_t).
  pair symbols:        A_{t^-1} = A_t^-1,  d_{t^-1} = A_t^-1 d_t        (T_{t^-1} = A_t^-1 T_t)
  involution symbols:  the relation must be symmetric: bar(A_s)^2 = id on V / T_s and A_s T_s = T_s.
     We store an exact involution A_s^2 = I with A_s d_s = d_s (every involution of V / T_s lifts to one; the
     lift graph only depends on bar(A_s)), so that A_{s^-1} = A_s^-1 = A_s and d_{s^-1} = d_s hold literally
     and every formula written for pairs applies verbatim with inv[s] = s.
L_t := T_{t^-1} = span(d_{t^-1}) is the 'backtrack line' of t at its source (walk t t^-1 adds L_t).

Chart-only conditions for a controller of girth >= 5 (every closed walk of length <= 4 backtracks):
  W4: reduced words of length 4 (t_{i+1} != t_i^-1; for an involution s this forbids 's s') controllable
      (columns A_{t4}A_{t3}A_{t2} d_{t1}, A_{t4}A_{t3} d_{t2}, A_{t4} d_{t3}, d_{t4} span F_p^m);
  C0: union over t != u of L_t + L_u and over (t, u), u != t^-1, of L_t + A_t^-1 L_u equals F_p^m
      (closed walks t t^-1 u u^-1 and t u u^-1 t^-1 from (v, x) cover x + those planes);
  C2: for each reduced 2-word (t, u): W = span(A_u d_t, d_u) (geodesic image), spur directions
      d_{r^-1} (spur at the end), A_u d_{r^-1} (spur at the middle vertex), A_u A_t d_{r^-1} (spur at the start)
      must hit every point of P(F^m / W) so the pencil of hyperplanes through W covers the fibre.
"""
import itertools, random, json
import rcl

class Alphabet:
    def __init__(self, p, m, inv, A, d):
        self.p, self.m = p, m
        self.K = len(inv)
        self.inv = list(inv)
        self.A = [[list(c) for c in a] for a in A]
        self.d = [list(v) for v in d]
        self.pairs = [t for t in range(self.K) if self.inv[t] > t]
        self.selfinv = [t for t in range(self.K) if self.inv[t] == t]

    # ---- consistency --------------------------------------------------------------------------------
    def check(self):
        p, m = self.p, self.m
        I = rcl.identity(m)
        for t in range(self.K):
            assert rcl.rank(self.A[t], p) == m, ('singular', t)
            assert any(self.d[t]), ('zero d', t)
        for t in self.pairs:
            ti = self.inv[t]
            assert rcl.matmul(self.A[ti], self.A[t], p) == I, ('inverse', t)
            assert rcl.matvec(self.A[ti], self.d[t], p) == self.d[ti], ('inverse d', t)
        for s in self.selfinv:
            assert rcl.matmul(self.A[s], self.A[s], p) == I, ('involution', s)
            assert rcl.matvec(self.A[s], self.d[s], p) == self.d[s], ('involution d', s)
        return True

    def set_pair(self, t, A, d):
        """set chart of pair symbol t and derive t^-1."""
        p = self.p; ti = self.inv[t]
        Ai = rcl.matinv(A, p)
        self.A[t] = [list(c) for c in A]; self.d[t] = list(d)
        self.A[ti] = Ai; self.d[ti] = rcl.matvec(Ai, d, p)

    # ---- words ---------------------------------------------------------------------------------------
    def reduced_words(self, L):
        for w in itertools.product(range(self.K), repeat=L):
            if all(w[i + 1] != self.inv[w[i]] for i in range(L - 1)): yield w

    def ctrl_columns(self, word):
        cols = []
        L = len(word)
        for i in range(L):
            v = list(self.d[word[i]])
            for j in range(i + 1, L): v = rcl.matvec(self.A[word[j]], v, self.p)
            cols.append(v)
        return cols

    def controllable(self, word):
        return rcl.rank(self.ctrl_columns(word), self.p) == self.m

    def bad_words(self, L=4):
        return [w for w in self.reduced_words(L) if not self.controllable(w)]

    # ---- chart-only conditions -----------------------------------------------------------------------
    def Lvec(self, t):
        """direction of L_t = T_{t^-1}."""
        return self.d[self.inv[t]]

    def c0_uncovered(self):
        p, m, K = self.p, self.m, self.K
        covered = set()
        for t in range(K):
            Lt = self.Lvec(t)
            Ati = self.A[self.inv[t]]
            for u in range(K):
                Lu = self.Lvec(u)
                covered |= {tuple(v) for v in rcl.span([Lt, Lu], p, m)}
                if u != self.inv[t]:
                    covered |= {tuple(v) for v in rcl.span([Lt, rcl.matvec(Ati, Lu, p)], p, m)}
        return p ** m - len(covered)

    def c2_missing(self, detail=False):
        p, m, K = self.p, self.m, self.K
        total = 0; det = {}
        for t in range(K):
            for u in range(K):
                if u == self.inv[t]: continue
                Au, At = self.A[u], self.A[t]
                W = [rcl.matvec(Au, self.d[t], p), self.d[u]]
                dimW = rcl.rank(W, p)
                Wset = rcl.span(W, p, m)
                dirs = []
                for r in range(K):
                    Lr = self.Lvec(r)
                    dirs.append(Lr)
                    v = rcl.matvec(Au, Lr, p); dirs.append(v)
                    dirs.append(rcl.matvec(Au, rcl.matvec(At, Lr, p), p))
                pts = set()
                for v in dirs:
                    if rcl.rank(W + [v], p) == dimW: continue
                    # canonical projective point of F^m / W: min over scalar multiples and coset
                    best = None
                    for l in range(1, p):
                        lv = [(l * x) % p for x in v]
                        for w in Wset:
                            c = tuple((a + b) % p for a, b in zip(lv, w))
                            if best is None or c < best: best = c
                    pts.add(best)
                nq = (p ** (m - dimW) - 1) // (p - 1)
                total += nq - len(pts)
                det[(t, u)] = (dimW, len(pts), nq)
        if detail: return total, det
        return total

    def score(self):
        return (len(self.bad_words(4)), self.c0_uncovered(), self.c2_missing())

    # ---- io -------------------------------------------------------------------------------------------
    def to_json(self):
        return {'p': self.p, 'm': self.m, 'inv': self.inv, 'A': self.A, 'd': self.d}

    @staticmethod
    def from_json(d):
        return Alphabet(d['p'], d['m'], d['inv'], d['A'], d['d'])

    def tsets(self):
        return [sorted(rcl.undigits(v, self.p) for v in rcl.span([self.d[t]], self.p, self.m)) for t in range(self.K)]

    def write_rcl(self, fn, N, out_edges, D=4):
        """rclift/oddlift input: out_edges[c] = list of (target, symbol)."""
        p, m = self.p, self.m; F = p ** m
        with open(fn, 'w') as f:
            f.write('%d %d %d %d %d\n' % (N, F, p, m, D))
            for c in range(N):
                f.write('%d' % len(out_edges[c]) + ''.join(' %d %d' % (dd, k) for dd, k in out_edges[c]) + '\n')
            f.write('%d\n' % self.K)
            ts = self.tsets()
            for k in range(self.K):
                f.write(' '.join(map(str, rcl.lin_table(self.A[k], p, m))) + '\n')
                f.write('%d ' % len(ts[k]) + ' '.join(map(str, ts[k])) + '\n')

# ---- random charts -------------------------------------------------------------------------------------
def random_vec(p, m, rng):
    while True:
        v = [rng.randrange(p) for _ in range(m)]
        if any(v): return v

def stabiliser_element(p, m, d, rng, small=False):
    """random P in GL_m(F_p) with P d = d (small: identity on all but one vector of a basis through d)."""
    # basis through d: d plus standard vectors completing it
    basis = [list(d)]
    for j in range(m):
        e = [0] * m; e[j] = 1
        if rcl.rank(basis + [e], p) == len(basis) + 1: basis.append(e)
    while True:
        imgs = [list(d)] + [list(b) for b in basis[1:]]
        if small:
            j = rng.randrange(1, m); imgs[j] = random_vec(p, m, rng)
        else:
            for j in range(1, m): imgs[j] = random_vec(p, m, rng)
        # P maps basis[i] -> imgs[i]; columns of P = P e_j: solve via matrices: P = IMG * BASIS^-1
        Binv = rcl.matinv(basis, p)   # columns of BASIS are basis vectors
        P = rcl.matmul(imgs, Binv, p)  # P e_j = IMG (BASIS^-1 e_j)
        if rcl.rank(P, p) == m: return P

def random_involution(p, m, d, rng):
    """random exact involution A with A d = d (over F_p; p = 2 handled via I + N, N^2 = 0, N d = 0)."""
    I = rcl.identity(m)
    if p != 2:
        # diagonal signs in a basis through d (d gets +1), random conjugate within the stabiliser of d
        basis = [list(d)]
        for j in range(m):
            e = [0] * m; e[j] = 1
            if rcl.rank(basis + [e], p) == len(basis) + 1: basis.append(e)
        signs = [1] + [rng.choice([1, p - 1]) for _ in range(m - 1)]
        if all(s == 1 for s in signs) and rng.random() < 0.9:
            signs[rng.randrange(1, m)] = p - 1
        imgs = [[(s * x) % p for x in b] for s, b in zip(signs, basis)]
        A0 = rcl.matmul(imgs, rcl.matinv(basis, p), p)
    else:
        while True:
            N = [[rng.randrange(2) for _ in range(m)] for _ in range(m)]
            if not any(any(c) for c in N): continue
            if rcl.matvec(N, d, 2) != [0] * m: continue
            if rcl.matmul(N, N, 2) != [[0] * m for _ in range(m)]: continue
            break
        A0 = [[(I[j][i] + N[j][i]) % 2 for i in range(m)] for j in range(m)]
    P = stabiliser_element(p, m, d, rng)
    Pinv = rcl.matinv(P, p)
    A = rcl.matmul(P, rcl.matmul(A0, Pinv, p), p)
    return A

def random_alphabet(p, m, npairs, nself, rng):
    K = 2 * npairs + nself
    inv = [0] * K
    for i in range(npairs): inv[2 * i] = 2 * i + 1; inv[2 * i + 1] = 2 * i
    for i in range(nself): inv[2 * npairs + i] = 2 * npairs + i
    A = [None] * K; d = [None] * K
    al = Alphabet(p, m, inv, [rcl.identity(m)] * K, [[1] + [0] * (m - 1)] * K)
    for i in range(npairs):
        al.set_pair(2 * i, rcl.random_gl(p, m, rng), random_vec(p, m, rng))
    for i in range(nself):
        s = 2 * npairs + i
        dv = random_vec(p, m, rng)
        al.d[s] = dv; al.A[s] = random_involution(p, m, dv, rng)
    return al

def load_extracted_chart(fn):
    """extract_chart.py style json (A columns modulo the translation line, Tb basis) -> exact Alphabet (repair columns)."""
    dd = json.load(open(fn))
    p, m = dd['p'], dd['m']
    inv = [dd['inv'][str(k)] for k in range(len(dd['A']))]
    A = []; d = []
    for k in range(len(dd['A'])):
        cols = [list(c) for c in dd['A'][k]]
        dv = list(dd['Tb'][k][0])
        # repair: add multiples of dv to columns until invertible
        for attempt in range(200):
            if rcl.rank(cols, p) == m: break
            j = attempt % m; l = 1 + (attempt // m) % (p - 1)
            cols[j] = [(a + l * b) % p for a, b in zip(cols[j], dv)]
        A.append(cols); d.append(dv)
    al = Alphabet(p, m, inv, A, d)
    # re-derive inverse symbols exactly from the lower index of each pair (the repaired columns may not be exact inverses)
    for t in al.pairs: al.set_pair(t, al.A[t], al.d[t])
    return al

if __name__ == '__main__':
    import sys
    dd = json.load(open(sys.argv[1]))
    al = load_extracted_chart(sys.argv[1]) if 'Tb' in dd else Alphabet.from_json(dd)
    al.check()
    print('symbols', al.K, 'inv', al.inv)
    print('score (bad words, C0 uncovered, C2 missing):', al.score())
    print('C2 detail:', al.c2_missing(detail=True)[1])

def word_index(al, L=4):
    """lexicographic index of reduced words (same order as ctrl5.c build_words)."""
    return {w: i for i, w in enumerate(al.reduced_words(L))}

def write_badset(al, fn, mode='w'):
    idx = word_index(al)
    bw = al.bad_words(4)
    with open(fn, mode) as f:
        f.write('%d %s\n' % (len(bw), ' '.join(str(idx[w]) for w in bw)))
    return len(bw)

def write_cht(al, fn):
    """chart tables in rclift/ctrl5 format: F p m K, then per symbol lin[F] and 't tau...'."""
    p, m = al.p, al.m; F = p ** m
    ts = al.tsets()
    with open(fn, 'w') as f:
        f.write('%d %d %d %d\n' % (F, p, m, al.K))
        for k in range(al.K):
            f.write(' '.join(map(str, rcl.lin_table(al.A[k], p, m))) + '\n')
            f.write('%d ' % len(ts[k]) + ' '.join(map(str, ts[k])) + '\n')
