#!/bin/bash
# 2026-11-28: integration test, escape_keep_rules + correlated decoder (EKR_CORR), d=4 1e-3 priors, small sample (own jsonl, not the record)
cd ~/Desktop/QC/prime_msd
export PYTHONPATH=.:scratchpad OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export EKR_CORR=1 EKR_PRIORS=runs/surface/priors_d4_s0.npz EKR_SEED=900000 EKR_SLICE=10000
export EKR_OUT=runs/surface/ekr_corr_test.jsonl EKR_CKPT=runs/surface/ckpt_corr_test EKR_CORR_DUMP=runs/surface/ekr_corr_test_s900000.npz
ulimit -v 3000000
nice -n 19 python3 -u scratchpad/escape_keep_rules.py 1e-3 400000 40000 5 par-hooks-f1 one+gap-tick1-closed
