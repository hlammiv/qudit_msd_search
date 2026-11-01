#!/bin/bash
# Sequential, single-thread candidate builds (p = 11, m = 4, r = 13) for the law validation,
# ~2 min each: 8 seeds at k = 2294 (26 checks) and 8 at k = 2296 (24 checks).
cd /home/hlamm/Desktop/QC/prime_msd/.claude/worktrees/agent-afd1ed4434f0e089f || exit 1
ulimit -v 4000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONPATH=.
for k in 2294 2296; do
  for s in 1 2 3 4 5 6 7 8; do
    f=scratchpad/routes/isd/runs/cands/p11_r13_k${k}_seed${s}.json
    [ -f "$f" ] || python3 scratchpad/rm_m4_candidate.py 11 13 $k $s "$f"
  done
done
echo ALLDONE
