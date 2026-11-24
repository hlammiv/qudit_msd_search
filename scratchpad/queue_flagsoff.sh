#!/bin/bash
# 2026-11-23: flags-off control for branch_t_cult (cat-flag detectors recorded but not post-selected): does the proxy gap appear?
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 PYTHONPATH=.:scratchpad
ulimit -v 3000000
o=runs/cultivation/branch_t
nice -n 19 python3 -u scratchpad/branch_t_cult.py run 2 1e-3 10000000 $o/flagsoff_color3_d2_p1e-3_s301.json 301 20000 par-hooks-f1 color3 flag > $o/flagsoff_color3_d2_p1e-3_s301.log 2>&1 &
nice -n 19 python3 -u scratchpad/branch_t_cult.py run 4 1e-3 10000000 $o/flagsoff_grow_d4_p1e-3_s401.json 401 20000 par-hooks-f1 grow flag > $o/flagsoff_grow_d4_p1e-3_s401.log 2>&1 &
nice -n 19 python3 -u scratchpad/branch_t_cult.py run 4 1e-3 10000000 $o/flagsoff_grow_d4_p1e-3_s402.json 402 20000 par-hooks-f1 grow flag > $o/flagsoff_grow_d4_p1e-3_s402.log 2>&1 &
wait
