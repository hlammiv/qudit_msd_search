#!/bin/bash
# 2026-11-03: f=1 escape variants end to end at 1e-3 (none vs one colour postselection, joint cuts), then 7e-4;
# alongside, the tick-idle logical CSUM at D = 9 (f = 1).  usage: queue_lenore_1103.sh
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1
PY=python3; [ -x ~/cap_validate_work/venv/bin/python ] && PY=~/cap_validate_work/venv/bin/python
echo "start $(date)"
(ulimit -v 4000000; nice -n 15 $PY -u -m primepower_msd.logical_csum sweep 4 4 9 1e-3,1.5e-3 400000 400000 seq tick1 >> runs/surface/csum_tick1_d9.log 2>&1) &
(ulimit -v 4000000; nice -n 15 $PY -u -m primepower_msd.logical_csum sweep 4 4 9 1e-3,1.5e-3 400000 400000 base tick1 >> runs/surface/csum_tick1_d9.log 2>&1) &
ULIMIT_KB=10000000 scratchpad/queue_e2e_p.sh 1e-3 16000000 10 par-hooks-f1-noneesc,par-hooks-f1-esc
ULIMIT_KB=10000000 scratchpad/queue_e2e_p.sh 7e-4 10000000 10 par-hooks,par-hooks-f1-esc
wait
echo "done $(date)"
