#!/bin/bash
# 2026-11-21 step 2: (A) d=8 priors 1e6 clean seed 200000 at 1e-3; (B) p=1.5e-3 d=4: seed-0 frames -> priors -> priors/no-priors runs (seed 100000).
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 PYTHONPATH=.:scratchpad
PY=python3; [ -x ~/cap_validate_work/venv/bin/python ] && PY=~/cap_validate_work/venv/bin/python
ulimit -v 4000000
V=one+gap-tick1-closed
case $1 in
A) EKR_D=8 EKR_SEED=200000 EKR_PRIORS=runs/surface/priors_d8_s0.npz nice -n 19 $PY -u scratchpad/escape_keep_rules.py 1e-3 4000000 1000000 12 par-hooks-f1 $V >> runs/surface/step2_d8pri.log 2>&1 ;;
B) F=runs/surface/frames_d4_p1.5e-3_s0.npz
   [ -f $F ] || EKR_FRAMES_ONLY=$F nice -n 19 $PY -u scratchpad/escape_keep_rules.py 1.5e-3 4000000 0 8 par-hooks-f1 $V >> runs/surface/step2_p15.log 2>&1
   ACC=$($PY -c "import numpy as np; print(int(np.load('$F')['accepted']))")
   [ -f runs/surface/priors_d4_p1.5e-3_s0.npz ] || nice -n 19 $PY -I scratchpad/make_priors.py $F $ACC 1.5e-3 1000000 runs/surface/priors_d4_p1.5e-3_s0.npz >> runs/surface/step2_p15.log 2>&1
   EKR_SEED=100000 EKR_PRIORS=runs/surface/priors_d4_p1.5e-3_s0.npz nice -n 19 $PY -u scratchpad/escape_keep_rules.py 1.5e-3 4000000 1000000 8 par-hooks-f1 $V >> runs/surface/step2_p15.log 2>&1
   EKR_SEED=100000 nice -n 19 $PY -u scratchpad/escape_keep_rules.py 1.5e-3 4000000 1000000 8 par-hooks-f1 $V >> runs/surface/step2_p15.log 2>&1 ;;
esac
echo "done $1 $(date)" >> runs/surface/step2_$1.done
