#!/bin/bash
# 2026-11-13 lenore: d = 16 residual-frame priors.  (1) seed-0 residual frames -> priors_d16_s0.npz; (2) closed window at
# seed 100000 (independent), without and with priors, concurrently, 6 workers each, 10k-shot slices (checkpointed).
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 PYTHONPATH=.:scratchpad
PY=~/cap_validate_work/venv/bin/python
ulimit -v 8000000
V=one+gap-tick1-closed
echo "frames $(date)"
EKR_D=16 EKR_FRAMES_ONLY=runs/surface/frames_d16_s0.npz nice -n 15 $PY -u scratchpad/escape_keep_rules.py 1e-3 4000000 1000 12 par-hooks-f1 $V >> runs/surface/escape_immediate_d16pri.log 2>&1
A=$($PY -I -c "import numpy as np; print(int(np.load('runs/surface/frames_d16_s0.npz')['accepted']))")
$PY -I scratchpad/make_priors.py runs/surface/frames_d16_s0.npz $A 1e-3 1000000 runs/surface/priors_d16_s0.npz >> runs/surface/escape_immediate_d16pri.log 2>&1
echo "validate $(date)"
EKR_LOGTAG=d16nopri EKR_SLICE=10000 EKR_SEED=100000 EKR_D=16 EKR_VARIANT=$V ULIMIT_KB=8000000 bash scratchpad/queue_immediate.sh 6 "1e-3 4000000 100000" &
EKR_LOGTAG=d16pri EKR_SLICE=10000 EKR_SEED=100000 EKR_D=16 EKR_VARIANT=$V ULIMIT_KB=8000000 EKR_PRIORS=runs/surface/priors_d16_s0.npz bash scratchpad/queue_immediate.sh 6 "1e-3 4000000 100000" &
wait
echo "done $(date)"
