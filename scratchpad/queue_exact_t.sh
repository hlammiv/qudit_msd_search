#!/bin/bash
# 2026-11-22: exact-T vs Clifford-proxy validation of stage-1 cultivation (ColorCult D=3), scratchpad/exact_t_cult.py.
# usage: queue_exact_t.sh "d p shots seed" ...   (one background process per spec; outputs runs/cultivation/exact_t/)
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 PYTHONPATH=.
ulimit -v 3000000
for spec in "$@"; do
  set -- $spec
  nice -n 19 python3 -u scratchpad/exact_t_cult.py run $1 $2 $3 runs/cultivation/exact_t/${CFG:-serial}_d$1_p$2_s$4.json $4 20000 ${CFG:-serial} > runs/cultivation/exact_t/${CFG:-serial}_d$1_p$2_s$4.log 2>&1 &
done
wait
