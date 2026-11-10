#!/bin/bash
# 2026-11-09: firm up the d = 4 + priors cut-7 point (1e-12 in one round rested on one event): new seed, 10x the clean sample.
cd "$(dirname "$0")/.."
EKR_SEED=200000 EKR_D=4 EKR_VARIANT=one+gap-tick1-closed EKR_PRIORS=runs/surface/priors_d4_s0.npz bash scratchpad/queue_immediate.sh 3 "1e-3 4000000 1000000"
