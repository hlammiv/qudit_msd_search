#!/bin/bash
# local 2026-11-02: d = 4 cultivation growth 3 -> 5 (r1c2), parallel schedule, idle noise at every tick f = 1 and 0.1;
# 1e7 shots per point at p = 5e-4, 1e-3; 2 workers, 20k-shot batches (~150 MB each; local memory is tight) -> runs/cultivation/sweep_grow35_d4_par_tick{1,0.1}.json
cd /home/hlamm/Desktop/QC/prime_msd
ulimit -v 1500000
export PYTHONPATH=. OMP_NUM_THREADS=1
echo "start $(date)"
for f in 1 0.1; do
  nice -n 19 python -u -m primepower_msd.cult_sweep --dims 4 --mode grow --rounds 1 --checks 2 --shots 10000000 --chunks 30 --jobs 2 --batch 20000 --ps 0.0005,0.001 --tick $f --parallel --out runs/cultivation/sweep_grow35_d4_par_tick$f.json >> runs/cultivation/cult_tick.log 2>&1
done
echo "done $(date)"
