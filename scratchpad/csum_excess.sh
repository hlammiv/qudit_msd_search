#!/bin/bash
# Direct measurement of the logical CSUM's excess over idling at fixed D (MEAS_DISTILL 5j): 1e6 shots per point,
# D = 5 (p = 1e-3, 2e-3) and D = 7 (p = 2e-3, 3e-3), both bases, Pauli-frame CSUM then the no-gate two-patch memory.
cd /home/hlamm/Desktop/QC/prime_msd
ulimit -v 4000000
export OMP_NUM_THREADS=1
for dec in seq base; do
  nice -n 19 python -u -m primepower_msd.logical_csum sweep 4 4 5 1e-3,2e-3 1000000 1000000 $dec >> runs/surface/csum_excess.log 2>&1
  nice -n 19 python -u -m primepower_msd.logical_csum sweep 4 4 7 2e-3,3e-3 1000000 1000000 $dec >> runs/surface/csum_excess.log 2>&1
done
echo "excess runs done $(date)"
