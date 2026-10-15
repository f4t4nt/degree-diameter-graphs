#!/bin/sh
# Positive control: rebuild the recorded voltage-lift stall (defect 16, i.e. 8 unordered pairs at
# distance 5), close it with the symmetry-breaking 2-opt repair, and verify the result by plain BFS.
# Usage: sh run_control.sh [REPAIR_SEED]   (default seed 1; every seed we tried reaches 0 in seconds)
set -e
SEED=${1:-1}
cc -O2 -o lift lift.c
cc -O2 -o repair repair.c -lm
cc -O2 -o verify_bfs ../../../tools/verify/verify_bfs.c
# SECONDS = 0: load the checkpoint and dump its graph as z8_b13_defect16.ck.edges
./lift cyc:8 13 0 29 0 z8_b13_defect16.ck
./verify_bfs z8_b13_defect16.ck.edges 4
rm -f repaired_104.edges
./repair z8_b13_defect16.ck.edges "$SEED" 60 repaired_104.edges
./verify_bfs repaired_104.edges 4
