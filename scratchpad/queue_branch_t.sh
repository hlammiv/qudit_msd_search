#!/bin/bash
# 2026-11-23: exact-T (branch-sum) vs proxy for d1=5 growth cultivation (GrowCult 3->5, par-hooks-f1), scratchpad/branch_t_cult.py.
# usage: queue_branch_t.sh "d p shots seed" ...   (one background process per spec; outputs runs/cultivation/branch_t/)
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 PYTHONPATH=.:scratchpad
ulimit -v 3000000
for spec in "$@"; do
  set -- $spec
  o=runs/cultivation/branch_t/grow_d$1_p$2_s$4
  nice -n 19 python3 -u scratchpad/branch_t_cult.py run $1 $2 $3 $o.json $4 20000 par-hooks-f1 grow > $o.log 2>&1 &
done
wait
