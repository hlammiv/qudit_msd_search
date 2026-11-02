#!/bin/bash
# local (lucia) 2026-11-01: d = 32 surface-code host sweep, same grid and shots as d = 4/8/16, 4 workers, nice 19,
# one line per finished point in runs/surface/zd_d32.jsonl.  Paused by polite_guard.py while others use > 12 cores.
cd /home/hlamm/Desktop/QC/prime_msd
ulimit -v 4000000
export OMP_NUM_THREADS=1
echo "start $(date)"
nice -n 19 python -u -m primepower_msd.surface_sweep 32 4 3,5,7,9 1e-3,1.5e-3,2e-3,3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,1e-2 100000 40000 >> runs/surface/zd_d32_sweep.log 2>&1
echo "done $(date)"
