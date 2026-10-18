#!/usr/bin/env python3
"""Per-edge route-chart lift over F_3^3 of a 5-regular controller (the family of R. Mizuno's (15,3) construction), stage 2 = type SAT
plus matrix realization plus lift (RESULTS.md section 2.4).

Chart family: at every vertex v the five departure lines D(v, i) (i = slot = position in adj[v]) are points of
PG(2,3): four of them on one projective line m_v ("ordinary" slots) and one off it ("special" slot).  A directed
edge e = (v1 -> v2) carries M_e in GL(3,3) with M_e D(e) = D(e^-1) (its departure line at v1 maps to the arrival
line at v2 = the departure line of the reverse dart) and the lift is (v1, x) ~ (v2, M_e x + F_3 d_e); the reverse
dart uses M_e^-1.  A 3-walk h e g (in-dart h at v1, out-dart g at v2) is full rank iff D(g) is not on the plane
M_e(span(D(h^-1), D(e))); with the coplanar configuration the bad walks through e are determined by its "type":
  ordinary-ordinary edge: M_e can send the two planes m_{v1}, span(sp_{v1}, D(e)) to the two planes through
      D(e^-1) that avoid the other departure points of v2: then NO walk through e (either direction) is bad and the
      edge pencils of e are complete;
  special (v1) - ordinary (v2): the 4 in-darts at v1 map bijectively to the 4 planes through D(e^-1) at v2:
      the one sent to m_{v2} has 3 bad continuations (the ordinary out-darts), the one sent to span(sp_{v2}, D(e^-1))
      has 1 (the special out-dart), the other two none; reverse direction: the 3 ordinary in-darts at v2 all have the
      bad continuation h_m^-1, the special in-dart has h_s^-1 (variables Tm, Ts: which in-dart plays which role);
  special-special: a bijection P between the in-darts at v1 and the out-darts at v2, bad iff P(h) = g (both ways).
SAT: choose the special slot per vertex and the types so that every ordered pair (u, w) of non-adjacent vertices and
every vertex (closed walk) keeps a good NB 3-walk; adjacent pairs are covered by the edge pencils.  A model is then
realized: coordinates per vertex (m_v = <e1, e2>, ordinary departure points = the 4 points of that line, special =
<e3>), M_e by brute force over GL(3,3) matching the required pencil map, the lift written as an edge list and checked
by a plain BFS here (and by the verifier in tools/ outside).
usage: typesat3.py <controller .ctrl|.edges> [--out lift.edges] [--time S] [--enum K] [--fix-special SLOT]
"""
import collections
import itertools
import sys
import time

from pysat.formula import IDPool
from pysat.solvers import Cadical153

P = 3


def load_controller(fn):
    adj = collections.defaultdict(list)
    if fn.endswith('.ctrl'):
        lines = [l for l in open(fn) if l.strip() and not l.startswith('#')]
        N, K = map(int, lines[0].split())
        for v in range(N):
            adj[v] = list(map(int, lines[3 + v].split()[:K]))
    else:
        for l in open(fn):
            if l.strip() and not l.startswith('#'):
                u, v = map(int, l.split()[:2]); adj[u].append(v); adj[v].append(u)
        N = 1 + max(adj)
    adj = [adj[v] for v in range(N)]
    for v in range(N):
        assert len(adj[v]) == 5 and len(set(adj[v])) == 5 and v not in adj[v], 'controller must be simple 5-regular'
    rev = [[adj[adj[v][i]].index(v) for i in range(5)] for v in range(N)]   # slot of the reverse dart
    return N, adj, rev


def walks3(N, adj, rev):
    """NB 3-walks as dart triples ((v0,i0),(v1,i1),(v2,i2)) grouped by (start, end)."""
    W = collections.defaultdict(list)
    for v0 in range(N):
        for i0 in range(5):
            v1 = adj[v0][i0]
            for i1 in range(5):
                if i1 == rev[v0][i0]: continue
                v2 = adj[v1][i1]
                for i2 in range(5):
                    if i2 == rev[v1][i1]: continue
                    v3 = adj[v2][i2]
                    W[(v0, v3)].append(((v0, i0), (v1, i1), (v2, i2)))
    return W


class TypeSAT:
    def __init__(self, N, adj, rev, fix_special=None):
        self.N, self.adj, self.rev = N, adj, rev
        self.pool = IDPool(); self.cnf = []
        pid = self.pool.id
        # special slot per vertex
        self.S = {(v, i): pid(('S', v, i)) for v in range(N) for i in range(5)}
        for v in range(N):
            self.cnf.append([self.S[(v, i)] for i in range(5)])
            for i in range(5):
                for j in range(i + 1, 5):
                    self.cnf.append([-self.S[(v, i)], -self.S[(v, j)]])
            if fix_special is not None:
                self.cnf.append([self.S[(v, fix_special)]])
        # per undirected edge (v1 < v2, slots i1 at v1, i2 at v2): type variables
        self.edges = {}
        for v1 in range(N):
            for i1 in range(5):
                v2 = adj[v1][i1]
                if v2 < v1: continue
                i2 = rev[v1][i1]
                ins1 = [h for h in range(5) if h != i1]   # darts at v1 other than e (their reverses are the in-darts)
                ins2 = [g for g in range(5) if g != i2]
                E = {'v1': v1, 'i1': i1, 'v2': v2, 'i2': i2, 'ins1': ins1, 'ins2': ins2}
                # case special at v1, ordinary at v2: Tm1[h], Ts1[h] over h in ins1 (exactly one each, distinct)
                E['Tm1'] = {h: pid(('Tm1', v1, i1, h)) for h in ins1}
                E['Ts1'] = {h: pid(('Ts1', v1, i1, h)) for h in ins1}
                # case ordinary at v1, special at v2: Tm2[g], Ts2[g] over g in ins2
                E['Tm2'] = {g: pid(('Tm2', v1, i1, g)) for g in ins2}
                E['Ts2'] = {g: pid(('Ts2', v1, i1, g)) for g in ins2}
                # case special-special: bijection Pm[h][g]
                E['P'] = {(h, g): pid(('P', v1, i1, h, g)) for h in ins1 for g in ins2}
                for T in ('Tm1', 'Ts1'):
                    self.exactly_one([E[T][h] for h in ins1])
                for T in ('Tm2', 'Ts2'):
                    self.exactly_one([E[T][g] for g in ins2])
                for h in ins1:
                    self.cnf.append([-E['Tm1'][h], -E['Ts1'][h]])
                for g in ins2:
                    self.cnf.append([-E['Tm2'][g], -E['Ts2'][g]])
                for h in ins1:
                    self.exactly_one([E['P'][(h, g)] for g in ins2])
                for g in ins2:
                    self.exactly_one([E['P'][(h, g)] for h in ins1])
                self.edges[(v1, i1)] = E
                self.edges[(v2, i2)] = E
        # bad-walk variables with the forcing implications
        self.B = {}
        self.W = walks3(N, adj, rev)
        for key, ws in self.W.items():
            for w in ws:
                self.bad_var(w)
        # coverage: non-adjacent ordered pairs and closed walks
        self.ncov = 0
        for u in range(N):
            for w in range(N):
                if w in adj[u]: continue
                ws = self.W.get((u, w), [])
                self.cnf.append([-self.B[x] for x in ws])
                self.ncov += 1

    def exactly_one(self, lits):
        self.cnf.append(list(lits))
        for a, b in itertools.combinations(lits, 2):
            self.cnf.append([-a, -b])

    def bad_var(self, w):
        if w in self.B: return self.B[w]
        (v0, i0), (v1, i1), (v2, i2) = w
        b = self.pool.id(('B',) + tuple(x for d in w for x in d)); self.B[w] = b
        E = self.edges[(v1, i1)]
        # orientation of the stored edge: stored v1 -> v2 with slots i1 (at E.v1) and i2 (at E.v2)
        forward = (E['v1'] == v1 and E['i1'] == i1)
        hslot = self.rev[v0][i0]        # slot at v1 of the reverse of the in-dart
        s1 = self.S[(v1, i1)]            # e's slot special at v1
        s2 = self.S[(v2, E['i2'] if forward else E['i1'])]   # reverse slot special at v2
        sh = self.S[(v1, hslot)]
        sg = self.S[(v2, i2)]
        if forward:
            Tm_here, Ts_here, Tm_there, Ts_there = E['Tm1'][hslot], E['Ts1'][hslot], E['Tm2'][i2], E['Ts2'][i2]
            Pv = E['P'][(hslot, i2)]
        else:
            Tm_here, Ts_here, Tm_there, Ts_there = E['Tm2'][hslot], E['Ts2'][hslot], E['Tm1'][i2], E['Ts1'][i2]
            Pv = E['P'][(i2, hslot)]
        # special here (v1), ordinary there (v2): in-dart h -> m_{v2} kills ordinary g; h -> sp-plane kills special g
        self.cnf.append([-s1, s2, -Tm_here, sg, b])
        self.cnf.append([-s1, s2, -Ts_here, -sg, b])
        # ordinary here, special there: ordinary in-darts killed at the out-dart g with Tm_there; special in-dart at Ts_there
        self.cnf.append([s1, -s2, sh, -Tm_there, b])
        self.cnf.append([s1, -s2, -sh, -Ts_there, b])
        # special-special: bijection
        self.cnf.append([-s1, -s2, -Pv, b])
        return b


# ---------- realization over F_3^3 ----------
def mv(M, x): return tuple(sum(M[i][j] * x[j] for j in range(3)) % 3 for i in range(3))


def all_gl3():
    """GL(3,3) as (M, image map point -> point)."""
    out = []
    pts = list(itertools.product(range(3), repeat=3))
    for cells in itertools.product(range(3), repeat=9):
        M = (cells[0:3], cells[3:6], cells[6:9])
        det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])) % 3
        if det: out.append((M, {x: mv(M, x) for x in pts}))
    return out


def mat_inv3(M):
    """inverse over F_3 by brute force on the adjugate: solve M X = I column by column."""
    pts = list(itertools.product(range(3), repeat=3))
    cols = []
    for j in range(3):
        e = tuple(int(i == j) for i in range(3))
        x = [p for p in pts if mv(M, p) == e][0]
        cols.append(x)
    return tuple(tuple(cols[j][i] for j in range(3)) for i in range(3))
def vadd(a, b): return tuple((x + y) % 3 for x, y in zip(a, b))
def smul(s, a): return tuple((s * x) % 3 for x in a)
def line(d): return frozenset(smul(s, d) for s in range(3))
def plane(a, b): return frozenset(vadd(smul(s, a), smul(t, b)) for s in range(3) for t in range(3))
def proj(v):
    for s in (1, 2):
        w = smul(s, v)
        if w > v: v = w
    return v


def realize(N, adj, rev, special, types, gl3):
    """special[v] = slot; types[(v1,i1)] for the stored orientation: ('oo',), ('so', hm, hs) (v1 special: in-dart slot
    hm -> m_{v2}, hs -> sp-plane), ('os', gm, gs) (v2 special), ('ss', perm dict h->g).  Returns M[(v,i)], d[(v,i)]."""
    e1, e2, e3 = (1, 0, 0), (0, 1, 0), (0, 0, 1)
    mpts = [e1, e2, (1, 1, 0), (1, 2, 0)]            # the four points of the line m = <e1, e2>
    D = {}
    for v in range(N):
        k = 0
        for i in range(5):
            if i == special[v]: D[(v, i)] = e3
            else: D[(v, i)] = mpts[k]; k += 1
    Mmap, dvec = {}, {}
    for (v1, i1), T in types.items():
        v2, i2 = adj[v1][i1], rev[v1][i1]
        De, Ae = D[(v1, i1)], D[(v2, i2)]
        sp1, sp2 = special[v1], special[v2]
        # required plane images: list of (plane at v1 through De, allowed set of planes at v2 through Ae)
        planes1 = {}
        for h in range(5):
            if h == i1: continue
            planes1[h] = plane(D[(v1, h)], De)
        planes2 = {}
        for g in range(5):
            if g == i2: continue
            planes2[g] = plane(D[(v2, g)], Ae)
        lineAe = line(Ae)
        ord1 = [h for h in range(5) if h != i1 and h != sp1]
        ord2 = [g for g in range(5) if g != i2 and g != sp2]
        m1 = plane(D[(v1, ord1[0])], De) if ord1 else None      # m_{v1} (through De) when De is ordinary
        s1p = plane(D[(v1, sp1)], De) if sp1 != i1 else None     # plane through sp_{v1} and De
        m2 = plane(D[(v2, ord2[0])], Ae) if ord2 else None
        s2p = plane(D[(v2, sp2)], Ae) if sp2 != i2 else None
        def ok(M, im):
            if line(im[De]) != lineAe: return False
            if T[0] == 'oo':
                i1m = frozenset(im[x] for x in m1); i1s = frozenset(im[x] for x in s1p)
                return i1m != m2 and i1m != s2p and i1s != m2 and i1s != s2p and i1m != i1s
            if T[0] == 'so':
                hm, hs = T[1], T[2]
                return frozenset(im[x] for x in planes1[hm]) == m2 and frozenset(im[x] for x in planes1[hs]) == s2p
            if T[0] == 'os':
                gm, gs = T[1], T[2]
                return frozenset(im[x] for x in m1) == planes2[gm] and frozenset(im[x] for x in s1p) == planes2[gs]
            if T[0] == 'ss':
                return all(frozenset(im[x] for x in planes1[h]) == planes2[g] for h, g in T[1].items())
            raise ValueError
        found = None
        for M, im in gl3:
            if ok(M, im): found = M; break
        if found is None: return None, (v1, i1, T)
        Mmap[(v1, i1)] = found
        dvec[(v1, i1)] = Ae
        Minv = mat_inv3(found)          # reverse dart: M^-1, d = M^-1 Ae (a generator of De)
        Mmap[(v2, i2)] = Minv
        dvec[(v2, i2)] = mv(Minv, Ae)
    return (Mmap, dvec), None


def build_lift(N, adj, Mmap, dvec):
    pts = list(itertools.product(range(3), repeat=3)); idx = {p: i for i, p in enumerate(pts)}
    edges = set()
    for v in range(N):
        for i in range(5):
            w = adj[v][i]; M = Mmap[(v, i)]; d = dvec[(v, i)]
            for x in pts:
                y0 = mv(M, x)
                for lam in range(3):
                    y = vadd(y0, smul(lam, d))
                    a, b = v * 27 + idx[x], w * 27 + idx[y]
                    if a == b: return None, 'loop'
                    edges.add((min(a, b), max(a, b)))
    return sorted(edges), None


def bfs_check(n, edges, D=3):
    adj = [[] for _ in range(n)]
    for a, b in edges: adj[a].append(b); adj[b].append(a)
    degs = sorted(set(len(x) for x in adj))
    far = 0
    for s in range(n):
        dist = [-1] * n; dist[s] = 0; fr = [s]; dd = 0
        while fr:
            dd += 1; nf = []
            for x in fr:
                for y in adj[x]:
                    if dist[y] < 0: dist[y] = dd; nf.append(y)
            fr = nf
        far += sum(1 for t in range(s + 1, n) if dist[t] < 0 or dist[t] > D)
    return degs, far


def main():
    args = sys.argv[1:]
    fn = args[0]
    def opt(name, default):
        return args[args.index(name) + 1] if name in args else default
    out = opt('--out', None); tlimit = float(opt('--time', '600')); enum = int(opt('--enum', '1'))
    fix = opt('--fix-special', None); fix = int(fix) if fix is not None else None
    t0 = time.time()
    N, adj, rev = load_controller(fn)
    enc = TypeSAT(N, adj, rev, fix_special=fix)
    print(f'TYPESAT N={N} vars={enc.pool.top} clauses={len(enc.cnf)} coverage_clauses={enc.ncov} build={time.time()-t0:.1f}s', flush=True)
    gl3 = None
    with Cadical153(bootstrap_with=enc.cnf) as s:
        found = 0
        while found < enum:
            s.conf_budget(int(2e6))
            r = s.solve_limited()
            if r is None:
                print(f'RESULT status=UNDECIDED time={time.time()-t0:.1f}s'); return
            if not r:
                print(f'RESULT status={"UNSAT" if found == 0 else "ENUM_COMPLETE"} models={found} time={time.time()-t0:.1f}s'); return
            found += 1
            model = set(l for l in s.get_model() if l > 0)
            special = {v: [i for i in range(5) if enc.S[(v, i)] in model][0] for v in range(N)}
            types = {}
            for key, E in enc.edges.items():
                if key != (E['v1'], E['i1']): continue
                v1, i1, v2, i2 = E['v1'], E['i1'], E['v2'], E['i2']
                s1, s2 = special[v1] == i1, special[v2] == i2
                if not s1 and not s2: types[(v1, i1)] = ('oo',)
                elif s1 and not s2:
                    hm = [h for h in E['ins1'] if E['Tm1'][h] in model][0]; hs = [h for h in E['ins1'] if E['Ts1'][h] in model][0]
                    types[(v1, i1)] = ('so', hm, hs)
                elif not s1 and s2:
                    gm = [g for g in E['ins2'] if E['Tm2'][g] in model][0]; gs = [g for g in E['ins2'] if E['Ts2'][g] in model][0]
                    types[(v1, i1)] = ('os', gm, gs)
                else:
                    perm = {h: [g for g in E['ins2'] if E['P'][(h, g)] in model][0] for h in E['ins1']}
                    types[(v1, i1)] = ('ss', perm)
            nsp = collections.Counter(t[0] for t in types.values())
            print(f'  model {found} at {time.time()-t0:.1f}s: special slots {collections.Counter(special.values())}, edge types {dict(nsp)}', flush=True)
            if gl3 is None: gl3 = all_gl3()
            real, err = realize(N, adj, rev, special, types, gl3)
            if real is None:
                print(f'  REALIZATION FAILED at {err}', flush=True)
                s.add_clause([-enc.S[(v, special[v])] for v in range(N)])
                continue
            Mmap, dvec = real
            edges, err = build_lift(N, adj, Mmap, dvec)
            if edges is None:
                print(f'  LIFT FAILED: {err}', flush=True); s.add_clause([-enc.S[(v, special[v])] for v in range(N)]); continue
            n = 27 * N
            degs, far = bfs_check(n, edges)
            print(f'  lift n={n} edges={len(edges)} degrees={degs} unordered pairs at distance > 3: {far} ({time.time()-t0:.1f}s)', flush=True)
            if out:
                with open(out, 'w') as f:
                    for a, b in edges: f.write(f'{a} {b}\n')
                print(f'  written {out}', flush=True)
                import json
                wit = {'controller_file': fn, 'N': N, 'adj': adj, 'special_slot': [special[v] for v in range(N)],
                       'types': {f'{v},{i}': (list(t) if t[0] != 'ss' else ['ss', {str(h): g for h, g in t[1].items()}]) for (v, i), t in types.items()},
                       'darts': {f'{v},{i}': {'to': adj[v][i], 'M': [list(r) for r in Mmap[(v, i)]], 'd': list(dvec[(v, i)])} for v in range(N) for i in range(5)},
                       'lift_rule': '(v, x) ~ (adj[v][i], M x + lambda d) for lambda in F_3; vertex index 27 v + 9 x0 + 3 x1 + x2 for x = (x0, x1, x2) in F_3^3',
                       'n': n, 'edges': len(edges), 'far_pairs': far}
                with open(out + '.witness.json', 'w') as f: json.dump(wit, f)
                print(f'  written {out}.witness.json', flush=True)
            if far == 0 and degs == [15]:
                print(f'LIFT_OK n={n} deg=15 D=3 file={out}', flush=True)
                return
            s.add_clause([-enc.S[(v, special[v])] for v in range(N)])
    print(f'RESULT status=DONE models={found} time={time.time()-t0:.1f}s')


if __name__ == '__main__':
    main()
