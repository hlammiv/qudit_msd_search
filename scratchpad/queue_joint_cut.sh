#!/bin/bash
# 2026-11-03: joint two-sector gap cut at f = 1 (escape_joint_cut.py), one+gap-tick1 vs none+gap-tick1, 3e5 trivial frames each
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 PYTHONPATH=.
ulimit -v 2600000
echo "start $(date)"
for v in none+gap-tick1 one+gap-tick1; do
  nice -n 15 python3 -u scratchpad/escape_joint_cut.py $v 300000 1e-3 3 >> runs/surface/escape_joint_cut.log 2>&1 &
done
wait
echo "done $(date)"
