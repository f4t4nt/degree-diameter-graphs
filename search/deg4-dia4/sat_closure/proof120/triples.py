"""Which family triples are UNSAT for every phi23 and every phi13 representative?  argv[1] = families."""
import sys
import time

from common import *

fams = sys.argv[1]
phis, auts, phi12, stab1, reps13 = setup()
t0 = time.time()
for rep, (phi13, orb) in enumerate(reps13):
    nsat = 0
    for phi23 in phis:
        if Reduced(phi12, phi13, phi23, families=fams).solve()[0]:
            nsat += 1
    print(f"families {fams}, phi13 rep {rep} (orbit {orb}): SAT for {nsat} of 6912 phi23", round(time.time() - t0, 1), "s", flush=True)
