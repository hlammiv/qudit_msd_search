"""Idle qudit-ticks of the d = 4 growth circuit (3 -> 5, r1c2) by phase and qudit role; fault-distance check on D = 3.
    PYTHONPATH=. python scratchpad/cult_idle_profile.py [schedule-kwargs as python dict]"""
import ast
import sys
from collections import defaultdict

from primepower_msd import cult_circuit as CC, cult_growth as CG
from primepower_msd.zdsim import add_tick_idle, check_determinism

kw = ast.literal_eval(sys.argv[1]) if len(sys.argv) > 1 else {"parallel": True}
marks = []


def wrap(cls, name):
    f = getattr(cls, name)

    def g(self, *a, **k):
        marks.append((len(self.c.ops), name + ("-final" if k.get("final") else "")))
        return f(self, *a, **k)
    setattr(cls, name, g)


for nm in ("encode", "stab_round", "check"):
    wrap(CC.ColorCult, nm)
wrap(CG.GrowCult, "grow")
g = CG.GrowCult(4, 1e-3, 1, 1, rounds3=1, checks3=2, rounds5=1, checks5=2, **kw)
nd = len(g.data5)
role = {q: ("data" if q < nd else "xa" if q < nd + 9 else "za" if q < nd + 18 else "cat") for q in range(g.c.nq)}
c2 = add_tick_idle(g.c, 1e-3)
base = set(id(o) for o in g.c.ops)
bounds = sorted(marks) + [(len(g.c.ops), "end")]
agg = defaultdict(float); k = seg = 0
for o in c2.ops:
    if id(o) in base:
        k += 1
        while seg + 1 < len(bounds) and k > bounds[seg + 1][0]:
            seg += 1
    elif o[0] == "DEP1":
        for q in o[1]:
            agg[(bounds[seg][1], role[q])] += o[2] / 1e-3
tot = defaultdict(float)
for (ph, r), v in agg.items():
    tot["data" if r == "data" else "anc"] += v
print("total idle", {k: round(v) for k, v in tot.items()}, "qudits", g.c.nq, "det", check_determinism(g.c, 64, 0)[:2])
print({k: round(v) for k, v in sorted(agg.items(), key=lambda kv: -kv[1]) if v >= 5})
cc = CC.ColorCult(3, 4, 1e-3, 1, 1, checks_per_round=2, rounds=1, **kw)
cc.c = add_tick_idle(cc.c, 1e-3)
print("D3 harmful single faults (tick idle):", cc.single_fault_scan()["undetected_harmful"])
