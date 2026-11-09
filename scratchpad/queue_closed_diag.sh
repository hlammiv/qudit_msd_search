#!/bin/bash
# 2026-11-09: closed-window gap diagnostics: per-shot residual dumps at d = 8 and d = 4 (same seeds as the 4e6 runs, few clean
# frames), then the d = 16 closed-window run at 1e-3.
cd "$(dirname "$0")/.."
V=one+gap-tick1-closed
EKR_DUMP=1 EKR_D=8 EKR_VARIANT=$V bash scratchpad/queue_immediate.sh 5 "1e-3 4000000 10000"
EKR_DUMP=1 EKR_D=4 EKR_VARIANT=$V bash scratchpad/queue_immediate.sh 5 "1e-3 4000000 10000"
EKR_D=16 EKR_VARIANT=$V bash scratchpad/queue_immediate.sh 5 "1e-3 4000000 1000000"
