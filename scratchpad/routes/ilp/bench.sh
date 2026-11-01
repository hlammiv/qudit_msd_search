#!/bin/bash
# Scaling table: solver vs MITM on random column subsets of A (p=7) and C (p=11).
# usage: bench.sh ENCODING TIME_LIMIT "A:100:5 A:200:5 ..." (instance:subset:w, subset 0 = full)
# Appends JSON lines to results.jsonl; one thread per run; 4 GB address-space cap.
cd "$(dirname "$0")/../../.." || exit 1
ulimit -v 4000000
export OMP_NUM_THREADS=1 PYTHONPATH=.
ENC=$1; TL=$2; shift 2
for spec in "$@"; do
  IFS=: read inst sub w <<< "$spec"
  if [ "$sub" = "0" ]; then subarg=""; else subarg="--subset $sub --seed 1"; fi
  python3 scratchpad/routes/ilp/ilp_distance.py --instance $inst $subarg --w $w --exact --encoding $ENC --time $TL \
      --log scratchpad/routes/ilp/results.jsonl --tag bench | cut -c1-160
done
