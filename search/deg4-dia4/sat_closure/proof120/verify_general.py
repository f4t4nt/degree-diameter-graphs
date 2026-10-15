"""Sharper statements behind the proof.

1. With kappa an ARBITRARY permutation of P1 (not necessarily an OA composite): (A)-(D) is UNSAT for
   every kappa != id tested (all 120 transpositions, all 3-cycles on a sample, 500 random permutations)
   and SAT for kappa = id.  Also (A)-(D) with kappa = id: print one solution.
2. Lemma 3 for the fixed phi12 and ALL 6912 phi13: phi13 o phi12^-1 is never an OA bijection
   (by the Aut x Aut symmetry this covers every pair of OA bijections).
3. Family subsets over all 6912 phi23 for rep 0 only for the pairs: counts of SAT instances, to make
   the remark on hypotheses precise (triples are covered by verify_proof.py's cores: ABCD cores of
   case (i) use all four families, case (ii) uses A, C, D).
"""
import itertools
import random
import sys
import time

from common import *

phis, auts, phi12, stab1, reps13 = setup()
phi13 = reps13[0][0]
phiset = set(phis)
inv12 = inverse(phi12)
t0 = time.time()

# 1. arbitrary kappa
ident = list(range(16))
m = Reduced(phi12, phi13, phis[0], kappa=ident)
r, model, _ = m.solve()
print("kappa = id: SAT =", r)
if r:
    a, rr, g = m.decode(model)
    print("  alpha =", a)
    print("  rho   =", rr)
    print("  sigma =", g)
random.seed(1)
tested = 0
for i, j in itertools.combinations(range(16), 2):
    k = list(ident)
    k[i], k[j] = k[j], k[i]
    assert not Reduced(phi12, phi13, phis[0], kappa=k).solve()[0], ("transposition SAT", i, j)
    tested += 1
for _ in range(500):
    k = list(ident)
    random.shuffle(k)
    if k == ident:
        continue
    assert not Reduced(phi12, phi13, phis[0], kappa=k).solve()[0], ("random kappa SAT", k)
    tested += 1
for i, j, l in itertools.combinations(range(16), 3):
    if random.random() < 0.1:
        k = list(ident)
        k[i], k[j], k[l] = k[j], k[l], k[i]
        assert not Reduced(phi12, phi13, phis[0], kappa=k).solve()[0], ("3-cycle SAT", i, j, l)
        tested += 1
print(f"arbitrary kappa != id: {tested} permutations tested, all UNSAT", round(time.time() - t0, 1), "s", flush=True)

# 2. Lemma 3 over all phi13
bad = 0
for phi13_ in phis:
    psi = tuple(phi13_[inv12[q]] for q in range(16))
    if psi in phiset:
        bad += 1
print("Lemma 3: phi13 o phi12^-1 is an OA bijection for", bad, "of 6912 phi13 (expect 0)")

# 3. pairs of families over all phi23, rep 0
for fams in ("AB", "AD", "BC", "BD", "CD"):
    nsat = 0
    for phi23 in phis:
        if Reduced(phi12, phi13, phi23, families=fams).solve()[0]:
            nsat += 1
    print(f"rep 0 families {fams}: SAT for {nsat} of 6912 phi23", round(time.time() - t0, 1), "s", flush=True)
