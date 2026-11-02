#!/bin/bash
# loranne 2026-11-02: logical CSUM (Pauli-frame decoder) vs no-gate two-patch memory with idle noise at every tick (f = 1),
# D = 5, 7 at p = 1e-3, 1.5e-3; 2e5 shots per point, 3 workers -> runs/surface/csum_d4_tick1.jsonl
cd ~/prime_msd_work
ulimit -v 2500000
export PYTHONPATH=pylib:. OMP_NUM_THREADS=1
echo "start $(date)"
for dec in seq base; do
  nice -n 15 python3 -u -m primepower_msd.logical_csum sweep 4 3 5,7 1e-3,1.5e-3 200000 200000 $dec tick1 >> runs/surface/csum_tick1.log 2>&1
done
echo "done $(date)"
