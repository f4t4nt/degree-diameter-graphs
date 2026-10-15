"""Check every step of the proof in PROOF120.md for one phi13 representative (argv[1] in 0,1,2) and
all 6912 phi23: for every cell x with z = kappa(x) != x, the exact constraint set named in the
corresponding case of the proof is unsatisfiable on its own.

  case (a) z on an OA line M of x:   {A[R(x);x], A[C(x);x], B[R(x);z], B[C(x);z], C[M;z]}
  case (b) z on a star S of x:       {C[L1(z);z], C[L2(z);z], D[L1(z);z], D[L2(z);z], A[S;x]}
  case (c) z on no line of x:        {A[R(x);x], A[C(x);x], B[R(x);z], B[C(x);z],
                                      C[L1(z);z], C[L2(z);z], D[L1(z);z], D[L2(z);z], A[R(z);x]}

Also: Lemma 3 (phi13 o phi12^-1 is not an OA bijection, so kappa != id for every phi23), and the
two transfer lemmas (Lemma 1: alpha(x) = kappa(x) from the four A/B constraints; Lemma 2:
rho^-1(x) = kappa(x) from the four C/D constraints) on the first 300 phi23.
"""
import sys
import time

from common import *

rep = int(sys.argv[1])
phis, auts, phi12, stab1, reps13 = setup()
phi13, orb = reps13[rep]
phiset = set(phis)
inv12 = inverse(phi12)
oa = [sorted(M) for M in oa_pms(phi12)]
row_of = {x: i for i, S in enumerate(ROWS) for x in S}
col_of = {x: 4 + i for i, S in enumerate(COLS) for x in S}
l1_of = {x: i for i, M in enumerate(oa[:4]) for x in M}
l2_of = {x: i for i, M in enumerate(oa[4:], start=4) for x in M}

psi = tuple(phi13[inv12[q]] for q in range(16))
assert psi not in phiset
print(f"rep {rep} (orbit {orb}): Lemma 3, phi13 o phi12^-1 is not an OA bijection: OK")

t0 = time.time()
n_case = {"a": 0, "b": 0, "c": 0}
n_id = 0
n_lemma = [0, 0]
for j, phi23 in enumerate(phis):
    m = Reduced(phi12, phi13, phi23)
    kappa = m.kappa
    if all(kappa[x] == x for x in range(16)):
        n_id += 1
        continue
    s = Cadical153(bootstrap_with=m.cl)
    for x in range(16):
        z = kappa[x]
        if z == x:
            continue
        AB = [("A", row_of[x], x), ("A", col_of[x], x), ("B", row_of[x], z), ("B", col_of[x], z)]
        CD = [("C", l1_of[z], z), ("C", l2_of[z], z), ("D", l1_of[z], z), ("D", l2_of[z], z)]
        on_star = row_of[x] == row_of[z] or col_of[x] == col_of[z]
        on_oa = l1_of[x] == l1_of[z] or l2_of[x] == l2_of[z]
        assert not (on_star and on_oa)  # fact F2: at most one of the four lines of x contains z
        if on_oa:
            M = l1_of[z] if l1_of[x] == l1_of[z] else l2_of[z]
            names = AB + [("C", M, z)]
            n_case["a"] += 1
        elif on_star:
            S = row_of[x] if row_of[x] == row_of[z] else col_of[x]
            names = CD + [("A", S, x)]
            n_case["b"] += 1
        else:
            names = AB + CD + [("A", row_of[z], x)]
            n_case["c"] += 1
        assert not s.solve(assumptions=[m.sel[n] for n in names]), (j, x, z, names)
        if j < 300:
            if not on_star:  # Lemma 1
                s2 = Cadical153(bootstrap_with=m.cl)
                s2.add_clause([-m.A[x][z]])
                assert not s2.solve(assumptions=[m.sel[n] for n in AB]), ("Lemma 1", j, x, z)
                s2.delete()
                n_lemma[0] += 1
            if not on_oa:  # Lemma 2
                s2 = Cadical153(bootstrap_with=m.cl)
                s2.add_clause([-m.R[z][x]])
                assert not s2.solve(assumptions=[m.sel[n] for n in CD]), ("Lemma 2", j, x, z)
                s2.delete()
                n_lemma[1] += 1
    s.delete()
    if j % 1000 == 999:
        print(j + 1, "phi23 done", round(time.time() - t0, 1), "s", flush=True)
print(f"rep {rep}: all 6912 phi23 checked; kappa=id instances: {n_id}; moved cells by case: {n_case}; "
      f"Lemma 1 / Lemma 2 checks: {n_lemma}; {round(time.time() - t0, 1)} s")
