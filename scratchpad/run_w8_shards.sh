#!/bin/bash
# usage: run_w8_shards.sh OUTDIR NW shard...   (runs the span-3 brute-force shards one after another)
out=$1; nw=$2; shift 2
for w in "$@"; do ${PY:-python3} scratchpad/rm20_w8_brute.py 3 $w $nw $out/span3_w$w.json > $out/span3_w$w.log 2>&1; done
