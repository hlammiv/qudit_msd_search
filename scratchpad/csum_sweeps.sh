#!/bin/bash
# logical CSUM sweeps, independent decoder (with baselines) then Pauli-frame decoder; 2 single-thread workers each
cd /home/hlamm/Desktop/QC/prime_msd
ulimit -v 4000000
export OMP_NUM_THREADS=1
nice -n 19 python -u -m primepower_msd.logical_csum sweep 4 2 3,5,7,9 1e-3,2e-3,3e-3,4e-3,6e-3 100000 40000 indep >> runs/surface/csum_sweep_indep.log 2>&1
nice -n 19 python -u -m primepower_msd.logical_csum sweep 4 2 3,5,7,9 1e-3,2e-3,3e-3,4e-3,6e-3 100000 40000 seq >> runs/surface/csum_sweep_seq.log 2>&1
echo "sweeps done $(date)"
