#!/bin/bash
# Sequential single-thread chain: Part 1 validation, exact restricted-column law test, Part 2.
cd /home/hlamm/Desktop/QC/prime_msd/.claude/worktrees/agent-afd1ed4434f0e089f || exit 1
ulimit -v 4000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export PYTHONPATH=.:scratchpad/routes/isd
R=scratchpad/routes/isd/runs
python3 -u scratchpad/routes/isd/run_validation.py $R p7a p7b p11
python3 -u scratchpad/routes/isd/law_restricted.py $R/p11_r13_k2290_seed1.npz $R/law_restricted_c13_s100.json 13 100 12 0
python3 -u scratchpad/routes/isd/law_restricted.py $R/p11_r13_k2290_seed1.npz $R/law_restricted_c14_s120.json 14 120 10 1
# wait for the candidate builder if it is still running
while ! grep -q ALLDONE $R/cands/build.log; do sleep 30; done
python3 -u scratchpad/routes/isd/law_part2.py $R/cands $R/law_part2.json 10 --exact4 2
echo CHAIN_DONE
