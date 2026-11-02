#!/bin/bash
# loranne 2026-11-01: d = 32 cultivation growth 3->5 (recipe of record r1 c2), 2e7 shots at p = 5e-4, 1e-3, 3 workers,
# chunks appended to the .partial.jsonl as they finish.  Paused by polite_guard.py while others use > 1 core.
cd ~/prime_msd_work
export OMP_NUM_THREADS=1 PYTHONPATH=pylib:.
ulimit -v 2500000
echo "start $(date)"
nice -n 15 python3 -u -m primepower_msd.cult_sweep --dims 32 --mode grow --rounds 1 --checks 2 --shots 20000000 --chunks 40 --jobs 3 --ps 0.0005,0.001 --out runs/cultivation/sweep_grow35_d32_long.json
echo "done $(date)"
