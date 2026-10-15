"""Independent second encoding of the 120 all-cross-double tightness problem (cross-check).

With every port holding one double of each of its two types, write phi_12, phi_13, phi_23 for the OA
bijections and transport everything to the ports of block 1:
  alpha: d12[p] ~ d13[alpha p]      (block-1 port of the {1,3}-neighbour)
  rho:   d12[p] ~ d23[phi12(rho p)] (the {2,3}-neighbour, named by the {1,2}-double at its 2-port)
  sigma: d13[p] ~ d23[phi12(sigma p)]
  kappa = phi13^-1 phi23 phi12 (fixed once phi23 is chosen).
Hub-to-port exactness for hubs of blocks 1 and 2 becomes, for every hub star S of block 1 and every
OA perfect matching M = phi12^-1(star of block 2):
  (A) S, rho S, alpha^-1 S, sigma S partition P1
  (B) S, alpha S, kappa rho S, kappa sigma S partition P1
  (C) M, alpha M, rho^-1 M, sigma^-1 M partition P1
  (D) kappa M, alpha M, kappa rho M, sigma^-1 M partition P1
(hubs of block 3 give two more families, not needed if these are already infeasible). The three
permutations are SAT variables; kappa is enumerated over all 6912 OA bijections phi23 for each of the
three phi13 representatives. This shares no code path with design120.Model120 beyond the OA list.
"""
import json
import sys
import time

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Cadical153

import design120 as D


def solve_reduced(phi12, phi13, phi23, families=("A", "B", "C", "D")):
    inv12 = {v: i for i, v in enumerate(phi12)}
    inv13 = {v: i for i, v in enumerate(phi13)}
    kappa = [inv13[phi23[phi12[p]]] for p in range(16)]
    kinv = {v: i for i, v in enumerate(kappa)}
    pool = IDPool()
    cl = []
    A = [[pool.id(("a", p, x)) for x in range(16)] for p in range(16)]
    R = [[pool.id(("r", p, x)) for x in range(16)] for p in range(16)]
    G = [[pool.id(("s", p, x)) for x in range(16)] for p in range(16)]
    T = pool.id(("T",))
    cl.append([T])

    def eq1(lits):
        for c in CardEnc.equals(lits=lits, bound=1, vpool=pool, encoding=EncType.seqcounter).clauses:
            cl.append(c)

    for M_ in (A, R, G):
        for p in range(16):
            eq1([M_[p][x] for x in range(16)])
            eq1([M_[x][p] for x in range(16)])
    stars = D.STARS
    oa_pms = [[inv12[q] for q in S] for S in stars]  # phi12^-1 of block-2 stars
    for S in stars:
        Sset = set(S)
        for x in range(16):
            if "A" in families:
                lits = ([T] if x in Sset else []) + [R[p][x] for p in S] + [A[x][p] for p in S] + [G[p][x] for p in S]
                eq1(lits)
            if "B" in families:
                y = kinv[x]
                lits = ([T] if x in Sset else []) + [A[p][x] for p in S] + [R[p][y] for p in S] + [G[p][y] for p in S]
                eq1(lits)
    for M in oa_pms:
        Mset = set(M)
        kM = {kappa[m] for m in M}
        for x in range(16):
            if "C" in families:
                lits = ([T] if x in Mset else []) + [A[m][x] for m in M] + [R[x][m] for m in M] + [G[x][m] for m in M]
                eq1(lits)
            if "D" in families:
                y = kinv[x]
                lits = ([T] if x in kM else []) + [A[m][x] for m in M] + [R[m][y] for m in M] + [G[x][m] for m in M]
                eq1(lits)
    s = Cadical153(bootstrap_with=cl)
    r = s.solve()
    s.delete()
    return r


def main():
    phis = D.oa_bijections()
    auts = D.automorphisms()
    phi12 = D.orbit_reps(phis, auts, auts)[0][0]
    stab1 = [a for a in auts if any(D.compose_phi(a, phi12, b) == phi12 for b in auts)]
    reps13 = D.orbit_reps(phis, stab1, auts)
    t0 = time.time()
    out = []
    for phi13, orb in reps13:
        counts = {"sat": 0, "unsat": 0}
        for phi23 in phis:
            if solve_reduced(phi12, phi13, phi23):
                counts["sat"] += 1
            else:
                counts["unsat"] += 1
        out.append({"phi13_orbit": orb, "phi23_tested": len(phis), **counts})
        print(out[-1], round(time.time() - t0, 1), "s", file=sys.stderr, flush=True)
    # sanity: single families must be satisfiable for at least one phi23
    sanity = {}
    for fam in (("A",), ("C",), ("A", "B"), ("C", "D")):
        sanity[",".join(fam)] = any(solve_reduced(phi12, reps13[0][0], phis[i], fam) for i in range(0, 6912, 691))
    print("sanity (families satisfiable alone):", sanity, file=sys.stderr, flush=True)
    json.dump({"results": out, "sanity": sanity}, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
