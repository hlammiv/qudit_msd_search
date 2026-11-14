#!/bin/bash
# 2026-11-14 lenore: d = 4 residual-frame priors at p = $1: seed-0 frames -> priors_d4_p$1_s0.npz -> closed window at seed
# 100000, without and with priors, concurrently, $2 workers each, 1e6 clean (checkpointed).
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 PYTHONPATH=.:scratchpad
PY=~/cap_validate_work/venv/bin/python
ulimit -v 4000000
V=one+gap-tick1-closed; P=$1; J=$2
echo "frames p=$P $(date)"
EKR_D=4 EKR_FRAMES_ONLY=runs/surface/frames_d4_p${P}_s0.npz nice -n 15 $PY -u scratchpad/escape_keep_rules.py $P 4000000 1000 $J par-hooks-f1 $V >> runs/surface/escape_immediate_pri_$P.log 2>&1
A=$($PY -I -c "import numpy as np; print(int(np.load('runs/surface/frames_d4_p${P}_s0.npz')['accepted']))")
$PY -I scratchpad/make_priors.py runs/surface/frames_d4_p${P}_s0.npz $A $P 1000000 runs/surface/priors_d4_p${P}_s0.npz >> runs/surface/escape_immediate_pri_$P.log 2>&1
echo "validate p=$P $(date)"
EKR_LOGTAG=nopri_$P EKR_SEED=100000 EKR_D=4 EKR_VARIANT=$V bash scratchpad/queue_immediate.sh $J "$P 4000000 1000000" &
EKR_LOGTAG=pri_$P EKR_SEED=100000 EKR_D=4 EKR_VARIANT=$V EKR_PRIORS=runs/surface/priors_d4_p${P}_s0.npz bash scratchpad/queue_immediate.sh $J "$P 4000000 1000000" &
wait
echo "done p=$P $(date)"
