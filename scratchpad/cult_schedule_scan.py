"""Scan stim-slot plaquette schedules for the parallel cultivation circuit: correctness (detector determinism), fault
distance (single-fault scan, D = 3 with tick idle) and idle qudit-ticks of the d = 4 growth circuit 3 -> 5 (r1c2).
    PYTHONPATH=. python scratchpad/cult_schedule_scan.py"""
import itertools as it
import sys

from primepower_msd.cult_circuit import ColorCult
from primepower_msd.cult_growth import GrowCult
from primepower_msd.zdsim import add_tick_idle, check_determinism


def idle(c, p=1e-3):
    c2 = add_tick_idle(c, p)
    return round((sum(len(o[1]) * o[2] for o in c2.ops if o[0] == "DEP1") - sum(len(o[1]) * o[2] for o in c.ops if o[0] == "DEP1")) / p)


cands = [("block-greedy", None)]
for ax, az, bz in it.product((1, 2), (1, 2), range(0, 13)):
    cands.append((f"X {ax}s, Z {az}s+{bz}", (ax, 0, az, bz)))
rows = []
for name, sch in cands:
    try:
        cc = ColorCult(3, 4, 1e-3, 1, 1, checks_per_round=2, rounds=1, parallel=True, schedule=sch)
        if check_determinism(cc.c, 64, 0):
            continue
        g = GrowCult(4, 1e-3, 1, 1, rounds3=1, checks3=2, rounds5=1, checks5=2, parallel=True, schedule=sch)
        if check_determinism(g.c, 64, 0):
            continue
        cc.c = add_tick_idle(cc.c, 1e-3)
        harm = cc.single_fault_scan()["undetected_harmful"]
        rows.append((idle(g.c), name, harm))
        print(name, "deterministic; D3 harmful (tick idle)", harm, "grow idle", rows[-1][0], flush=True)
    except Exception as e:
        print(name, "error", repr(e)[:80])
print("best:", sorted(r for r in rows if r[2] == 0)[:5])
