#!/bin/bash
# 2026-11-03: end-to-end cultivation -> escape at a given p.  usage: queue_e2e_p.sh P SHOTS JOBS CONFIGS
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1
[ -d pylib ] && export PYTHONPATH=pylib:. || export PYTHONPATH=.
PY=python3; [ -x ~/cap_validate_work/venv/bin/python ] && PY=~/cap_validate_work/venv/bin/python
ulimit -v 2600000
echo "start $(date)"
nice -n 15 $PY -u scratchpad/cult_escape_e2e.py $1 $2 $3 $4 >> runs/cultivation/cult_escape_e2e.log 2>&1
echo "done $(date)"
