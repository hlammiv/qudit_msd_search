#!/bin/bash
# Local queue: wait for the Y47 census, then run checkpointed censuses for t=45 and t=40, each under the guard.
cd /home/hlamm/Desktop/QC/prime_msd
while pgrep -f "[w]8_mitm44v2.py paired runs/weighted/grow_Y47" >/dev/null; do sleep 60; done
for t in 45 40; do
  setsid bash -c "ulimit -v 4000000; MITM_COORDS=4 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=12 PYTHONPATH=.:scratchpad exec nice -n 19 python -u scratchpad/w8_mitm44v2.py paired runs/weighted/pc_smoke.json $t 0 runs/weighted/mitm44/v2local_dpp_t$t.json" > runs/weighted/mitm44/v2local_dpp_t$t.log 2>&1 < /dev/null &
  sleep 5
  S=$(ps -eo sid,cmd | grep "[w]8_mitm44v2.py paired runs/weighted/pc_smoke.json $t" | grep -v bash | awk '{print $1}' | head -1)
  GUARD_MEM_MIN_GB=1.5 python3 -u scratchpad/polite_guard_mem.py $S 6 3 30 > runs/weighted/mitm44/v2local_guard_t$t.log 2>&1 < /dev/null
  echo "$(date) t=$t finished" >> runs/weighted/mitm44/queue_local.log
done
