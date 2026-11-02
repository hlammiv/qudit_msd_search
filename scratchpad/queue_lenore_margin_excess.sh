#!/bin/bash
# lenore 2026-11-01, two concurrent jobs:
# A. escape gap-cut margin: gap recipe at p = 1e-3, 6 more replicas x 4e6 attempts per sector, 12 workers
#    (per-task results -> graft_variants.jsonl; full gap curve + failing gaps per task -> graft_variants_gap_margin.log)
# B. p-dependence of the CSUM excess at fixed D: Pauli-frame CSUM vs no-gate two-patch memory, 1e6 shots, 8 workers
cd ~/cap_validate_work/prime_msd
PY=~/cap_validate_work/venv/bin/python
export PYTHONPATH=. OMP_NUM_THREADS=1
echo "start $(date)"
( nice -n 10 $PY -u -m primepower_msd.graft_run variants - 4000000 1e-3 12 one+pairs-cap+gap 1 4 6 >> runs/surface/graft_variants_gap_margin.log 2>&1; echo "margin done $(date)" ) &
for dec in seq base; do
  nice -n 10 $PY -u -m primepower_msd.logical_csum sweep 4 8 5 1.5e-3,3e-3,4e-3 1000000 1000000 $dec >> runs/surface/csum_excess_p.log 2>&1
  nice -n 10 $PY -u -m primepower_msd.logical_csum sweep 4 8 7 1.5e-3,4e-3 1000000 1000000 $dec >> runs/surface/csum_excess_p.log 2>&1
  nice -n 10 $PY -u -m primepower_msd.logical_csum sweep 4 8 9 3e-3,4e-3 1000000 1000000 $dec >> runs/surface/csum_excess_p.log 2>&1
  echo "excess $dec done $(date)"
done
wait
echo "all done $(date)"
