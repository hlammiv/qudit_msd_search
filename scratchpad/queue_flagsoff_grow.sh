#!/bin/bash
# 2026-11-24: flags-off control at d1=5 (GrowCult d=4, 1e-3) with the branch cap raised to 3M rows (was 400k: 18 overflows / 1.6e5)
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 PYTHONPATH=.:scratchpad BT_CAP=3000000
ulimit -v 7000000
o=runs/cultivation/branch_t
for s in 401 402 403; do
  nice -n 19 python3 -u scratchpad/branch_t_cult.py run 4 1e-3 10000000 $o/flagsoff_grow_d4_p1e-3_s$s.json $s 20000 par-hooks-f1 grow flag > $o/flagsoff_grow_d4_p1e-3_s$s.log 2>&1 &
done
wait
