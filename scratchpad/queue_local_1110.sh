#!/bin/bash
# 2026-11-10: local only (lenore full), one run at a time, 5 workers: (1) d = 16 closed-window re-run (OOM-killed 11-09 when
# stacked), (2) d = 4 + priors with 1e7 clean frames (settle the 1e-12 one-round point; clean floor).
cd "$(dirname "$0")/.."
V=one+gap-tick1-closed
EKR_D=16 EKR_VARIANT=$V bash scratchpad/queue_immediate.sh 5 "1e-3 4000000 1000000"
EKR_SEED=300000 EKR_D=4 EKR_VARIANT=$V EKR_PRIORS=runs/surface/priors_d4_s0.npz bash scratchpad/queue_immediate.sh 5 "1e-3 4000000 10000000"
