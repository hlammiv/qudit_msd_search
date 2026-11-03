#!/bin/bash
# wait for the escape tick floor queue (out file $1) to finish, then run the end-to-end configs $2 with $3 jobs under a guard (12 / 6)
cd "$(dirname "$0")/.."
until grep -q "^done" "$1" 2>/dev/null; do sleep 120; done
setsid nohup scratchpad/queue_e2e.sh $3 $2 > runs/cultivation/queue_e2e_esc.out 2>&1 < /dev/null &
sleep 3
S=$(ps -eo sid,args | grep "[q]ueue_e2e.sh $3 $2" | awk '{print $1}' | head -1)
echo "e2e SID $S"
[ -n "$S" ] && setsid nohup python3 scratchpad/polite_guard.py $S $4 $5 60 > runs/cultivation/guard_e2e_esc.log 2>&1 < /dev/null &
