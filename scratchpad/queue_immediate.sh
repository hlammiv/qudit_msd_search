#!/bin/bash
# 2026-11-07: immediate-decision escape (one colour class, cut 0) + syndrome-only gap, f = 1.
# env EKR_VARIANT (default one+gap-tick1-syn).
# usage: queue_immediate.sh JOBS "p shots n_clean" ["p shots n_clean" ...]
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 PYTHONPATH=.:scratchpad${EXTRA_PYPATH:+:$EXTRA_PYPATH}
PY=python3; [ -x ~/cap_validate_work/venv/bin/python ] && PY=~/cap_validate_work/venv/bin/python
ulimit -v ${ULIMIT_KB:-4000000}
J=$1; shift
for spec in "$@"; do
  set -- $spec
  echo "start $spec $(date)"
  nice -n 15 $PY -u scratchpad/escape_keep_rules.py $1 $2 $3 $J par-hooks-f1 ${EKR_VARIANT:-one+gap-tick1-syn} >> runs/surface/escape_immediate.log 2>&1
  echo "done $spec $(date)"
done
