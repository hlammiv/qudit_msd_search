#!/bin/bash
# 2026-11-02: escape (one+gap) with idle noise at every tick, p = 1e-3, both bases.
# usage: queue_escape_tick.sh JOBS SHOTS VARIANT REPLICAS ULIMIT_KB   (local: 3 500000 one+gap-tick1 3; loranne: 1 200000 one+gap-tick0.1 1)
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1
[ -d pylib ] && export PYTHONPATH=pylib:. || export PYTHONPATH=.
ulimit -v $5
echo "start $(date)"
nice -n 15 python3 -u -m primepower_msd.graft_run variants - $2 1e-3 $1 $3 1 4 $4 >> runs/surface/graft_variants_tick.log 2>&1
echo "done $(date)"
