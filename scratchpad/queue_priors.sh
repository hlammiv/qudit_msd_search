#!/bin/bash
# 2026-11-09: residual-frame priors in the escape decoder (priors from the seed-0 dumps), validated on an independent
# cultivation sample (EKR_SEED = 100000): closed window with and without priors on the same shots, d = 4 then d = 8.
cd "$(dirname "$0")/.."
V=one+gap-tick1-closed
for D in 4 8; do
  EKR_SEED=100000 EKR_D=$D EKR_VARIANT=$V bash scratchpad/queue_immediate.sh 3 "1e-3 4000000 100000"
  EKR_SEED=100000 EKR_D=$D EKR_VARIANT=$V EKR_PRIORS=runs/surface/priors_d${D}_s0.npz bash scratchpad/queue_immediate.sh 3 "1e-3 4000000 100000"
done
