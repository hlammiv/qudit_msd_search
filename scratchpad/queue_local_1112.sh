#!/bin/bash
# 2026-11-12: d = 4 + priors, seed 300000, 1e7 clean, moved from lenore (100/406 slices checkpointed there); 3 workers local.
cd "$(dirname "$0")/.."
EKR_LOGTAG=prior EKR_SEED=300000 EKR_D=4 EKR_VARIANT=one+gap-tick1-closed EKR_PRIORS=runs/surface/priors_d4_s0.npz bash scratchpad/queue_immediate.sh 3 "1e-3 4000000 10000000"
