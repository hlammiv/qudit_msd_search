#!/bin/bash
# local 2026-11-02: d = 16 escape, the d = 8 comparison set (variants one, one+gap; p = 1e-3, 2e-3; 2e5 attempts per sector and
# point), moved here from lenore while lenore is busy.  3 workers (~1.5 GB, ~6 h per task; 8 tasks).
cd /home/hlamm/Desktop/QC/prime_msd
ulimit -v 6000000
export PYTHONPATH=. OMP_NUM_THREADS=1
echo "start $(date)"
nice -n 19 python -u -m primepower_msd.graft_run variants - 200000 1e-3,2e-3 3 one,one+gap 1 16 >> runs/surface/graft_variants_d16.log 2>&1
echo "done $(date)"
