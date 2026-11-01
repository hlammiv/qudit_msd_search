#!/bin/bash
# Two concurrent lenore jobs at 6 workers each (2026-10-30):
# A. trimmed escape floor: gap recipe at p = 1e-3, 3 replicas x 4e6 attempts per sector (results per task -> graft_variants.jsonl)
# B. cultivation check-count sweep: growth 3->5, d = 4, (rounds, checks) = (1,1), (1,3), (2,2) vs the recipe of record (1,2),
#    1e8 shots per point at p = 5e-4 and 1e-3, one output file per variant
cd ~/cap_validate_work/prime_msd
PY=~/cap_validate_work/venv/bin/python
export PYTHONPATH=. OMP_NUM_THREADS=1
echo "start $(date)"
( nice -n 10 $PY -u -m primepower_msd.graft_run variants - 4000000 1e-3 6 one+pairs-cap+gap 1 4 3 >> runs/surface/graft_variants_gap_floor.log 2>&1; echo "floor done $(date)" ) &
for rc in "1 1" "1 3" "2 2"; do
  set -- $rc
  nice -n 10 $PY -u -m primepower_msd.cult_sweep --dims 4 --Ds 3 --mode grow --rounds $1 --checks $2 --shots 100000000 --chunks 50 --jobs 6 --ps 0.0005,0.001 --out runs/cultivation/sweep_grow35_d4_r$1c$2.json >> runs/cultivation/sweep_grow35_d4_checks.log 2>&1
  echo "checks r$1c$2 done $(date)"
done
wait
echo "all done $(date)"
