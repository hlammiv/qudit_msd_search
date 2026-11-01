#!/bin/bash
# d = 16 checks at 6 workers (2026-10-31): host-code memory sweep on the d = 8 grid, then growth cultivation 3->5 at 2e7 shots
cd ~/cap_validate_work/prime_msd
PY=~/cap_validate_work/venv/bin/python
export PYTHONPATH=. OMP_NUM_THREADS=1
echo "start $(date)"
nice -n 10 $PY -u -m primepower_msd.surface_sweep 16 8 3,5,7,9 1e-3,1.5e-3,2e-3,3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,1e-2 100000 40000 >> runs/surface/zd_d16_sweep.log 2>&1
echo "host d16 done $(date)"
nice -n 10 $PY -u -m primepower_msd.cult_sweep --dims 16 --Ds 3 --mode grow --rounds 1 --checks 2 --shots 20000000 --chunks 40 --jobs 8 --ps 0.0005,0.001 --out runs/cultivation/sweep_grow35_d16_long.json >> runs/cultivation/sweep_grow35_d16_long.log 2>&1
echo "growth d16 done $(date)"
