#!/bin/bash
# 2026-11-14 lenore: $1 = p7 (resume the 7e-4 with/without-priors pair from its checkpoints, 4 + 4) |
# d16 (d = 16 + priors and without, seed 200000, 1e6 clean, 3 + 3, 10k slices).
cd "$(dirname "$0")/.."
V=one+gap-tick1-closed
case $1 in
 p7)  EKR_LOGTAG=nopri_7e-4 EKR_SEED=100000 EKR_D=4 EKR_VARIANT=$V bash scratchpad/queue_immediate.sh 4 "7e-4 4000000 1000000" &
      EKR_LOGTAG=pri_7e-4 EKR_SEED=100000 EKR_D=4 EKR_VARIANT=$V EKR_PRIORS=runs/surface/priors_d4_p7e-4_s0.npz bash scratchpad/queue_immediate.sh 4 "7e-4 4000000 1000000" & ;;
 d16) EKR_LOGTAG=d16nopri_b EKR_SLICE=10000 EKR_SEED=200000 EKR_D=16 EKR_VARIANT=$V ULIMIT_KB=8000000 bash scratchpad/queue_immediate.sh 3 "1e-3 4000000 1000000" &
      EKR_LOGTAG=d16pri_b EKR_SLICE=10000 EKR_SEED=200000 EKR_D=16 EKR_VARIANT=$V ULIMIT_KB=8000000 EKR_PRIORS=runs/surface/priors_d16_s0.npz bash scratchpad/queue_immediate.sh 3 "1e-3 4000000 1000000" & ;;
esac
wait; echo "done $1 $(date)"
