"""Sweep the Z_d surface-code memory experiment over (D, basis, p); append JSON lines to runs/surface/zd_d{d}.jsonl.

python -m primepower_msd.surface_sweep d jobs [Dlist] [plist] [shots_low shots_high]
"""
from __future__ import annotations

import json
import sys
import time

from joblib import Parallel, delayed

from .surface_zd import Experiment


def point(D, d, basis, p, shots, seed):
    t = time.time()
    ex = Experiment(D, d, basis)
    r = ex.run(p, shots, seed)
    r["seconds"] = round(time.time() - t, 1)
    return r


def main():
    d = int(sys.argv[1]); jobs = int(sys.argv[2])
    Ds = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [3, 5, 7, 9]
    ps = [float(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else [1e-3, 1.5e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 8e-3, 1e-2]
    lo, hi = (int(sys.argv[5]), int(sys.argv[6])) if len(sys.argv) > 6 else (100_000, 40_000)
    tasks = [(D, basis, p, lo if p <= 2e-3 else hi) for D in Ds for basis in "ZX" for p in ps]
    tasks.sort(key=lambda t: -t[0])
    out = f"runs/surface/zd_d{d}.jsonl"
    res = Parallel(n_jobs=jobs, verbose=0)(delayed(point)(D, d, b, p, s, hash((D, b, p)) % 2**31) for D, b, p, s in tasks)
    with open(out, "a") as f:
        for r in res:
            f.write(json.dumps(r) + "\n")
    print("wrote", len(res), "points to", out)


if __name__ == "__main__":
    main()
