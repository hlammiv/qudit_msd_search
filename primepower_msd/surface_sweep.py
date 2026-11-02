"""Sweep the Z_d surface-code memory experiment over (D, basis, p); append JSON lines to runs/surface/zd_d{d}.jsonl.

python -m primepower_msd.surface_sweep d jobs [Dlist] [plist] [shots_low shots_high] [tick[f]]
(tick f: DEP1(f p) idle noise at every tick, surface_zd.memory_circuit(tick_idle=f), f = 1 by default -> runs/surface/zd_d{d}_tick{f}.jsonl)
"""
from __future__ import annotations

import json
import sys
import time

from joblib import Parallel, delayed

from .surface_zd import Experiment


def point(D, d, basis, p, shots, seed, tick=0.0):
    t = time.time()
    ex = Experiment(D, d, basis, tick_idle=tick)
    r = ex.run(p, shots, seed)
    r["seconds"] = round(time.time() - t, 1)
    return r


def main():
    d = int(sys.argv[1]); jobs = int(sys.argv[2])
    Ds = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [3, 5, 7, 9]
    ps = [float(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else [1e-3, 1.5e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 8e-3, 1e-2]
    lo, hi = (int(sys.argv[5]), int(sys.argv[6])) if len(sys.argv) > 6 else (100_000, 40_000)
    tick = float(sys.argv[7][4:] or 1) if len(sys.argv) > 7 and sys.argv[7].startswith("tick") else 0.0   # tick, tick0.1
    tasks = [(D, basis, p, lo if p <= 2e-3 else hi) for D in Ds for basis in "ZX" for p in ps]
    tasks.sort(key=lambda t: -t[0])
    out = f"runs/surface/zd_d{d}{f'_tick{tick:g}' if tick else ''}.jsonl"
    n = 0
    for r in Parallel(n_jobs=jobs, verbose=0, return_as="generator_unordered")(
            delayed(point)(D, d, b, p, s, hash((D, b, p)) % 2**31, tick) for D, b, p, s in tasks):
        with open(out, "a") as f:                  # one line per finished point: a crash or kill loses only the points in flight
            f.write(json.dumps(r) + "\n")
        n += 1
        print(r["D"], r["basis"], r["p"], "fails", r["fails"], "shots", r["shots"], flush=True)
    print("wrote", n, "points to", out)


if __name__ == "__main__":
    main()
