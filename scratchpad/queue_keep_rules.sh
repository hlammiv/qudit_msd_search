#!/bin/bash
# 2026-11-04: keep rules of the f = 1 escape (lighter vs recipe), 1e-3: residual frames of 1.6e7 cultivation shots + 4e6 clean frames.
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 PYTHONPATH=.:scratchpad
PY=python3; [ -x ~/cap_validate_work/venv/bin/python ] && PY=~/cap_validate_work/venv/bin/python
ulimit -v 6000000
echo "start $(date)"
nice -n 15 $PY -u scratchpad/escape_keep_rules.py 1e-3 16000000 4000000 20 par-hooks-f1 none+gap-tick1,one+gap-tick1 >> runs/surface/escape_keep_rules.log 2>&1
echo "done $(date)"
