"""Expected duration of a cultivation attempt and of an escape attempt with early abort (2026-11-08).
An attempt stops at the tick of its first firing postselected detector (cultivation: every detector; escape: the colour-region
detectors, rule none); a surviving attempt runs to the end (escape: end of the syndrome rounds, decision before the readout).
Ticks from the ASAP schedule of add_tick_idle (one tick per gate / CSUM pair).  python scratchpad/attempt_time.py shots"""
import json
import sys

import numpy as np

from cult_escape_e2e import parse
from primepower_msd.cult_growth import GrowCult
from primepower_msd.graft_run import VARIANTS
from primepower_msd.graft_zd import Escape
from primepower_msd.surface_zd import P_BUILD
from primepower_msd.zdsim import FrameSim, add_tick_idle


def det_ticks(c):
    t, mt = [0] * c.nq, []
    for op in c.ops:
        n = op[0]
        if n in ("F", "FD", "S", "R", "M"):
            for q in op[1]:
                t[q] += 1
                if n == "M":
                    mt.append(t[q])
        elif n == "CSUM":
            for a, b in op[1]:
                s = max(t[a], t[b]) + 1; t[a] = t[b] = s
    assert len(mt) == c.n_meas
    return np.array([max(mt[i] for i, _ in terms) if terms else 0 for terms in c.detectors]), max(t)


def abort_ticks(Dv, ticks, end):
    fire = Dv != 0
    first = np.where(fire, ticks[None, :], np.iinfo(np.int64).max).min(axis=1)
    return np.where(fire.any(axis=1), first, end)


shots = int(sys.argv[1])
kw, f, _ = parse("par-hooks-f1")
ex = Escape(4, "Z", **VARIANTS["one+gap-tick1"])
etk, edepth = det_ticks(ex.mem.circ)
post = ex.mem.det_post
tags = ex.mem.det_tags
t_end = max(etk[i] for i, tg in enumerate(tags) if not tg.startswith("final"))   # last syndrome-round detector
rows = []
for p in (1e-3, 7e-4, 5e-4):
    g = GrowCult(4, p, 1, 1, rounds3=1, rounds5=1, **kw)
    ctk, cdepth = det_ticks(g.c)
    c = add_tick_idle(g.c, f * p)
    s = FrameSim(c, shots, np.random.default_rng(11), gauge=False).run()
    Dv = s.detector_values()
    ta = abort_ticks(Dv, ctk, cdepth)
    acc = float((~np.any(Dv != 0, axis=1)).mean())
    se = FrameSim(ex.mem.circ, shots, np.random.default_rng(12)).run(scale=p / P_BUILD)
    De = se.detector_values()[:, post]
    te = abort_ticks(De, etk[post], t_end)
    keep_mask = float((~np.any(De != 0, axis=1)).mean())
    r = dict(p=p, shots=shots, cult_depth=int(cdepth), cult_accept=acc, cult_mean_ticks=float(ta.mean()),
             cult_mean_ticks_rejected=float(ta[ta < cdepth].mean()) if (ta < cdepth).any() else None,
             esc_depth=int(edepth), esc_decision_tick=int(t_end), esc_mask_keep=keep_mask, esc_mean_ticks=float(te.mean()))
    rows.append(r); print(json.dumps(r), flush=True)
with open("runs/surface/attempt_time.jsonl", "a") as fh:
    for r in rows:
        fh.write(json.dumps(r) + "\n")
