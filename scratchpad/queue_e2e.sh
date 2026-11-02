#!/bin/bash
# 2026-11-02: end-to-end cultivation -> escape (scratchpad/cult_escape_e2e.py) at p = 1e-3, 4e6 cultivation shots per config.
# usage: queue_e2e.sh JOBS CONFIGS   (loranne: serial,serial-hooks,par-hooks; local: par-hooks-f0.1,par-hooks-f1)
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1
[ -d pylib ] && export PYTHONPATH=pylib:. || export PYTHONPATH=.
ulimit -v 2500000
echo "start $(date)"
nice -n 15 python3 -u scratchpad/cult_escape_e2e.py 1e-3 4000000 $1 $2 >> runs/cultivation/cult_escape_e2e.log 2>&1
echo "done $(date)"
