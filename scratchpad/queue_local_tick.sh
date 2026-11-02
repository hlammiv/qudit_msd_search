#!/bin/bash
# local 2026-11-02: d = 4 surface memory with idle noise at every tick, strength f = 1 (SD6-like) and f = 0.1 (SI1000-like);
# 2 workers each -> runs/surface/zd_d4_tick{1,0.1}.jsonl
cd /home/hlamm/Desktop/QC/prime_msd
ulimit -v 3000000
export PYTHONPATH=. OMP_NUM_THREADS=1
echo "start $(date)"
nice -n 19 python -u -m primepower_msd.surface_sweep 4 2 3,5,7,9 5e-4,7e-4,1e-3,1.5e-3,2e-3,3e-3,4e-3,5e-3 100000 40000 tick1 >> runs/surface/zd_d4_tick1_sweep.log 2>&1 &
nice -n 19 python -u -m primepower_msd.surface_sweep 4 2 3,5,7,9 1e-3,1.5e-3,2e-3,3e-3,4e-3,5e-3,6e-3,8e-3 100000 40000 tick0.1 >> runs/surface/zd_d4_tick0.1_sweep.log 2>&1 &
wait
echo "done $(date)"
