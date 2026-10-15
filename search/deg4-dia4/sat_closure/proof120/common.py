"""Shared structure for the PROOF120 exploration: phi12 representative, the three phi13
representatives, the affine-plane classes behind phi12, and a selector-based encoding of the
reduced puzzle so that unsatisfiable cores can be extracted with assumptions.

design120.py lives one directory up (sat_closure/); it is put on the path here.
"""
import itertools
import pathlib
import sys

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Cadical153

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import design120 as D  # noqa: E402

PORTS = D.PORTS
PIDX = D.PIDX
STARS = D.STARS  # 0..3 rows, 4..7 columns
ROWS = STARS[:4]
COLS = STARS[4:]


def inverse(perm):
    inv = [0] * len(perm)
    for i, v in enumerate(perm):
        inv[v] = i
    return inv


def setup():
    """phi12 representative, stabiliser, phi13 representatives (as in design120_reduced.main)."""
    phis = D.oa_bijections()
    auts = D.automorphisms()
    phi12 = D.orbit_reps(phis, auts, auts)[0][0]
    stab1 = [a for a in auts if any(D.compose_phi(a, phi12, b) == phi12 for b in auts)]
    reps13 = D.orbit_reps(phis, stab1, auts)
    return phis, auts, phi12, stab1, reps13


def oa_pms(phi):
    """The 8 OA perfect matchings of phi (preimages of the 8 stars): 4 L1-classes then 4 L2-classes."""
    inv = inverse(phi)
    return [frozenset(inv[q] for q in S) for S in STARS]


def fifth_class(phi12):
    """The fifth parallel class: cells pairwise in different rows, columns, L1- and L2-classes."""
    lines = [frozenset(S) for S in STARS] + oa_pms(phi12)

    def collinear(x, y):
        return any(x in L and y in L for L in lines)

    classes = []
    seen = set()
    for x in range(16):
        if x in seen:
            continue
        cl = frozenset(y for y in range(16) if y == x or not collinear(x, y))
        classes.append(cl)
        seen |= cl
    return classes


class Reduced:
    """Reduced puzzle with one selector literal per (family, star-or-PM index, cell) constraint.

    Constraint names: (fam, i, x) where fam in "ABCD", i indexes the star (A, B) or the OA PM (C, D)
    and x the cell of P1 that must be covered exactly once.
    """

    def __init__(self, phi12, phi13, phi23, families="ABCD", hubfree=False, kappa=None):
        inv12 = inverse(phi12)
        if kappa is None:
            inv13 = inverse(phi13)
            kappa = [inv13[phi23[phi12[p]]] for p in range(16)]
        self.kappa = list(kappa)
        kinv = inverse(self.kappa)
        self.pool = pool = IDPool()
        self.cl = cl = []
        self.A = A = [[pool.id(("a", p, x)) for x in range(16)] for p in range(16)]
        self.R = R = [[pool.id(("r", p, x)) for x in range(16)] for p in range(16)]
        self.G = G = [[pool.id(("s", p, x)) for x in range(16)] for p in range(16)]
        T = pool.id(("T",))
        cl.append([T])
        self.sel = {}

        def eq1(lits, name=None):
            if name is None:
                for c in CardEnc.equals(lits=lits, bound=1, vpool=pool, encoding=EncType.seqcounter).clauses:
                    cl.append(c)
            else:
                s = pool.id(("sel",) + name)
                self.sel[name] = s
                for c in CardEnc.equals(lits=lits, bound=1, vpool=pool, encoding=EncType.seqcounter).clauses:
                    cl.append([-s] + c)

        for M_ in (A, R, G):
            for p in range(16):
                eq1([M_[p][x] for x in range(16)])
                eq1([M_[x][p] for x in range(16)])
        if hubfree:  # optional: alpha, rho, sigma map every cell off its row and column
            for S in STARS:
                for p in S:
                    for x in S:
                        for M_ in (A, R, G):
                            cl.append([-M_[p][x]])
        self.oa = oa = [sorted(M) for M in oa_pms(phi12)]
        for i, S in enumerate(STARS):
            Sset = set(S)
            for x in range(16):
                if "A" in families:
                    lits = ([T] if x in Sset else []) + [R[p][x] for p in S] + [A[x][p] for p in S] + [G[p][x] for p in S]
                    eq1(lits, ("A", i, x))
                if "B" in families:
                    y = kinv[x]
                    lits = ([T] if x in Sset else []) + [A[p][x] for p in S] + [R[p][y] for p in S] + [G[p][y] for p in S]
                    eq1(lits, ("B", i, x))
        for i, M in enumerate(oa):
            Mset = set(M)
            kM = {self.kappa[m] for m in M}
            for x in range(16):
                if "C" in families:
                    lits = ([T] if x in Mset else []) + [A[m][x] for m in M] + [R[x][m] for m in M] + [G[x][m] for m in M]
                    eq1(lits, ("C", i, x))
                if "D" in families:
                    y = kinv[x]
                    lits = ([T] if x in kM else []) + [A[m][x] for m in M] + [R[m][y] for m in M] + [G[x][m] for m in M]
                    eq1(lits, ("D", i, x))

    def solve(self, names=None, extra=()):
        """SAT status with the named constraints switched on (default: all)."""
        if names is None:
            names = list(self.sel)
        s = Cadical153(bootstrap_with=self.cl)
        for c in extra:
            s.add_clause(c)
        r = s.solve(assumptions=[self.sel[n] for n in names])
        core = None
        model = None
        if r:
            model = s.get_model()
        else:
            core = s.get_core()
        s.delete()
        return r, model, core

    def decode(self, model):
        pos = set(l for l in model if l > 0)
        a = [next(x for x in range(16) if self.A[p][x] in pos) for p in range(16)]
        r = [next(x for x in range(16) if self.R[p][x] in pos) for p in range(16)]
        g = [next(x for x in range(16) if self.G[p][x] in pos) for p in range(16)]
        return a, r, g

    def shrink(self, names):
        """Deletion-based minimal unsatisfiable subset of the given (unsat) constraint set."""
        names = list(names)
        s = Cadical153(bootstrap_with=self.cl)
        assert not s.solve(assumptions=[self.sel[n] for n in names])
        core = set(n for n in names if self.sel[n] in set(s.get_core()))
        names = [n for n in names if n in core]
        i = 0
        while i < len(names):
            trial = names[:i] + names[i + 1:]
            if s.solve(assumptions=[self.sel[n] for n in trial]):
                i += 1
            else:
                core = set(s.get_core())
                names = [n for n in trial if self.sel[n] in core]
        s.delete()
        return names


def cell(x):
    return PORTS[x]
