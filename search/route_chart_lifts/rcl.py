#!/usr/bin/env python3
"""Helpers for route-chart lifts: finite-field vectors as base-p integers, charts, rclift input files.

A chart = (A, T): A an m x m matrix over F_p given by its columns (A e_i), T a translation subspace
given by a basis (list of vectors).  A symbol alphabet has pairs (t, t^-1) with A_{t^-1} = A_t^-1 and
T_{t^-1} = A_t^-1 T_t, and self-inverse symbols with A^2 = I, A T = T.
"""
import itertools, random, subprocess, os

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

def matvec(cols, v, p):
    m = len(v); r = [0] * m
    for i in range(m):
        if v[i]:
            for j in range(m): r[j] = (r[j] + v[i] * cols[i][j]) % p
    return r

def matmul(A, B, p):
    """columns of A*B: A applied to each column of B."""
    return [matvec(A, b, p) for b in B]

def matinv(cols, p):
    m = len(cols)
    # rows of the matrix: M[j][i] = cols[i][j]; augment with identity, row reduce
    M = [[cols[i][j] % p for i in range(m)] + [1 if i == j else 0 for i in range(m)] for j in range(m)]
    for c in range(m):
        piv = next((r for r in range(c, m) if M[r][c]), None)
        if piv is None: return None
        M[c], M[piv] = M[piv], M[c]
        inv = pow(M[c][c], -1, p)
        M[c] = [(inv * u) % p for u in M[c]]
        for r in range(m):
            if r != c and M[r][c]:
                f = M[r][c]; M[r] = [(u - f * w) % p for u, w in zip(M[r], M[c])]
    # inverse matrix rows are M[j][m:]; return columns
    return [[M[j][m + i] for j in range(m)] for i in range(m)]

def rank(vecs, p):
    rows = [list(v) for v in vecs]
    if not rows: return 0
    m = len(rows[0]); r = 0
    for col in range(m):
        piv = next((i for i in range(r, len(rows)) if rows[i][col] % p), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][col], -1, p); rows[r] = [(inv * u) % p for u in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][col] % p:
                f = rows[i][col]; rows[i] = [(u - f * w) % p for u, w in zip(rows[i], rows[r])]
        r += 1
    return r

def span(basis, p, m):
    S = [[0] * m]
    for b in basis:
        S = [[(s[j] + l * b[j]) % p for j in range(m)] for s in S for l in range(p)]
    return S

def lin_table(cols, p, m):
    F = p ** m
    return [undigits(matvec(cols, digits(x, p, m), p), p) for x in range(F)]

def random_gl(p, m, rng):
    while True:
        cols = [[rng.randrange(p) for _ in range(m)] for _ in range(m)]
        if rank(cols, p) == m: return cols

def identity(m):
    return [[1 if i == j else 0 for j in range(m)] for i in range(m)]

class Alphabet:
    """symbols 0..K-1 with inv[k]; charts A[k] (columns), Tb[k] (basis of translation space)."""
    def __init__(self, p, m, pairs, selfinv, A, Tb):
        self.p, self.m = p, m
        self.K = 2 * pairs + selfinv
        self.inv = {}
        for i in range(pairs): self.inv[2 * i] = 2 * i + 1; self.inv[2 * i + 1] = 2 * i
        for i in range(selfinv): self.inv[2 * pairs + i] = 2 * pairs + i
        self.A, self.Tb = A, Tb

    def tsets(self):
        return [sorted(undigits(v, self.p) for v in span(self.Tb[k], self.p, self.m)) for k in range(self.K)]

    def controllability_ok(self, word):
        """word = list of symbols t1..tL (walk order); M = [A_tL..A_t2 b_t1 | ... | b_tL] columns from T's."""
        cols = []
        L = len(word)
        for i in range(L):
            vecs = [list(b) for b in self.Tb[word[i]]]
            for j in range(i + 1, L): vecs = [matvec(self.A[word[j]], v, self.p) for v in vecs]
            cols += vecs
        return rank(cols, self.p) == self.m

    def reduced_words(self, L):
        for w in itertools.product(range(self.K), repeat=L):
            if all(w[i + 1] != self.inv[w[i]] for i in range(L - 1)): yield w

    def count_bad_words(self, L, limit=None):
        bad = 0; total = 0
        for w in self.reduced_words(L):
            total += 1
            if not self.controllability_ok(w):
                bad += 1
                if limit is not None and bad >= limit: break
        return bad, total

def random_alphabet(p, m, pairs, selfinv, r, rng):
    """random chart satisfying the inverse conditions; self-inverse symbols get A = involution with A T = T."""
    A = [None] * (2 * pairs + selfinv); Tb = [None] * (2 * pairs + selfinv)
    for i in range(pairs):
        M = random_gl(p, m, rng)
        while True:
            B = [[rng.randrange(p) for _ in range(m)] for _ in range(r)]
            if rank(B, p) == r: break
        Minv = matinv(M, p)
        A[2 * i] = M; Tb[2 * i] = B
        A[2 * i + 1] = Minv; Tb[2 * i + 1] = [matvec(Minv, b, p) for b in B]
    for i in range(selfinv):
        k = 2 * pairs + i
        # involution: conjugate of a diagonal +-1 matrix (p odd) or a random involution built from a transvection-like form
        while True:
            P = random_gl(p, m, rng); Pinv = matinv(P, p)
            Dg = [[(rng.choice([1, p - 1]) if i_ == j else 0) for j in range(m)] for i_ in range(m)] if p != 2 else None
            if p == 2:
                # random involution over F_2: I + N with N^2 = 0; take N = u v^T with v^T u = 0
                u = [rng.randrange(2) for _ in range(m)]; v = [rng.randrange(2) for _ in range(m)]
                if not any(u) or not any(v) or sum(a * b for a, b in zip(u, v)) % 2: continue
                Dg = [[(1 if i_ == j else 0) ^ (u[i_] & v[j]) for i_ in range(m)] for j in range(m)]  # columns
            M = matmul(P, matmul(Dg, Pinv, p), p)
            # translation space: must be A-invariant: pick T inside an eigenspace-ish: use r random vectors v, replace by v + A v (A-invariant for involutions when p=2) or eigenvectors
            cand = []
            tries = 0
            while len(cand) < r and tries < 200:
                tries += 1
                v = [rng.randrange(p) for _ in range(m)]
                Av = matvec(M, v, p)
                w = [(a + b) % p for a, b in zip(v, Av)]  # eigenvector for +1 (A w = A v + v = w)
                if any(w) and rank(cand + [w], p) == len(cand) + 1: cand.append(w)
            if len(cand) == r: break
        A[k] = M; Tb[k] = cand
    return Alphabet(p, m, pairs, selfinv, A, Tb)

def write_input(fn, N, out_edges, alph, D=4):
    """out_edges[c] = list of (d, symbol)."""
    p, m = alph.p, alph.m; F = p ** m
    with open(fn, 'w') as f:
        f.write('%d %d %d %d %d\n' % (N, F, p, m, D))
        for c in range(N):
            f.write('%d' % len(out_edges[c]) + ''.join(' %d %d' % (d, k) for d, k in out_edges[c]) + '\n')
        f.write('%d\n' % alph.K)
        ts = alph.tsets()
        for k in range(alph.K):
            f.write(' '.join(map(str, lin_table(alph.A[k], p, m))) + '\n')
            f.write('%d ' % len(ts[k]) + ' '.join(map(str, ts[k])) + '\n')

RCLIFT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bin', 'rclift')

def run_rclift(fn, extra=()):
    out = subprocess.run([RCLIFT, fn, '--quiet'] + list(extra), capture_output=True, text=True, check=True).stdout.strip()
    d = {}
    for tok in out.split():
        k, v = tok.split('=')
        d[k] = v
    return d
