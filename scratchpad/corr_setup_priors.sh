#!/bin/bash
# 2026-11-28: pickled correlated-decoder setup with the d=4 1e-3 residual-frame priors (both bases), for escape_keep_rules EKR_CORR=1
cd ~/Desktop/QC/prime_msd
export PYTHONPATH=.:scratchpad OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ulimit -v 6000000
for b in Z X; do
  nice -n 19 python3 -u -c "import time,corr_escape as CE; t=time.time(); CE.cached_setup(4,'$b',1e-3,'one+gap-tick1-closed','runs/surface/priors_d4_s0.npz'); print('$b done', round(time.time()-t),'s',flush=True)"
done
