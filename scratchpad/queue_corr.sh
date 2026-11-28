#!/bin/bash
# 2026-11-27: correlated vs base escape decoder, clean frames, d=4, 1e-3, 1e6 shots per basis (checkpointed; rerun resumes)
cd ~/Desktop/QC/prime_msd
export PYTHONPATH=.:scratchpad OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ulimit -v 3000000
nice -n 19 python3 -u scratchpad/corr_escape.py 1e-3 1000000 6 3 4 runs/surface/corr_escape_d4_p1e-3_s3.npz
