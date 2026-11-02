#!/bin/bash
# lenore 2026-11-02: escape at d = 16 and 32, the d = 8 comparison set (variants one, one+gap; p = 1e-3, 2e-3; 2e5 attempts
# per sector and point; no pair certification).  8 workers each (d = 16 ~1.5 GB / 6 h per task, d = 32 ~5.1 GB / 26 h).
# Per-task results -> runs/surface/graft_variants_d{16,32}.jsonl; summaries + gap curves -> runs/surface/graft_variants_d{16,32}.log
cd ~/cap_validate_work/prime_msd
PY=~/cap_validate_work/venv/bin/python
export PYTHONPATH=. OMP_NUM_THREADS=1
echo "start $(date)"
for d in 16 32; do
  ( nice -n 10 $PY -u -m primepower_msd.graft_run variants - 200000 1e-3,2e-3 8 one,one+gap 1 $d >> runs/surface/graft_variants_d$d.log 2>&1; echo "d=$d done $(date)" ) &
done
wait
echo "all done $(date)"
