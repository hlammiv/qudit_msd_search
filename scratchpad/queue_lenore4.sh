#!/bin/bash
# Lenore queue at 4 workers (the box sits at load ~26/32 from other users; no load gate):
# 1. d = 8 escape, 2e5 attempts per sector, one-class recipe with and without the gap, p = 1e-3 and 2e-3
# 2. escape floor under the gap recipe at p = 1e-3: 25 replicas x 4e6 attempts per sector (1e8 per sector)
# 3. d = 8 logical CSUM sweep, Pauli-frame decoder (D = 3..9, five p)
# 4. d = 8 cultivation growth 3->5, 2e7 shots per point at p = 5e-4 and 1e-3
cd ~/cap_validate_work/prime_msd
PY=~/cap_validate_work/venv/bin/python
export PYTHONPATH=. OMP_NUM_THREADS=1
echo "start $(date)"
nice -n 10 $PY -u -m primepower_msd.graft_run variants - 200000 1e-3,2e-3 4 one,one+gap 1 8 1 >> runs/surface/graft_variants_d8.log 2>&1
echo "job1 d8 escape done $(date)"
nice -n 10 $PY -u -m primepower_msd.graft_run variants - 4000000 1e-3 4 one+pairs-cap+gap 1 4 25 >> runs/surface/graft_variants_gap_1e8.log 2>&1
echo "job2 1e8 floor done $(date)"
nice -n 10 $PY -u -m primepower_msd.logical_csum sweep 8 4 3,5,7,9 1e-3,2e-3,3e-3,4e-3,6e-3 100000 40000 seq >> runs/surface/csum_sweep_d8_seq.log 2>&1
echo "job3 d8 csum done $(date)"
nice -n 10 $PY -u -m primepower_msd.cult_sweep --dims 8 --Ds 3 --mode grow --rounds 1 --checks 2 --shots 20000000 --chunks 40 --jobs 4 --ps 0.0005,0.001 --out runs/cultivation/sweep_grow35_d8_long.json >> runs/cultivation/sweep_grow35_d8_long.log 2>&1
echo "job4 d8 cultivation done $(date)"
