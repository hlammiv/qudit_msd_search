#!/bin/bash
# 2026-11-11 lenore (free): $1 = d16 | prior.  Per-slice checkpoints (runs/surface/ckpt), progress lines in the log.
cd "$(dirname "$0")/.."
V=one+gap-tick1-closed
case $1 in
  d16)   ULIMIT_KB=8000000 EKR_D=16 EKR_VARIANT=$V bash scratchpad/queue_immediate.sh 8 "1e-3 4000000 1000000" ;;
  prior) EKR_SEED=300000 EKR_D=4 EKR_VARIANT=$V EKR_PRIORS=runs/surface/priors_d4_s0.npz bash scratchpad/queue_immediate.sh 10 "1e-3 4000000 10000000" ;;
esac
