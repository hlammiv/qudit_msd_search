"""Lenore sweep: stage-1 colour-code cultivation (Clifford proxy) for d in {2,4,8}, D in {3,5}, several p.

    PYTHONPATH=. python -m primepower_msd.cult_sweep --jobs 30 --shots 20000000 --out runs/cultivation/sweep.json
"""
from __future__ import annotations

import argparse
import json
import os
import time

from .cult_circuit import ColorCult


def one(d, D, p, shots, seed, rounds, checks):
    cc = ColorCult(D, d, p, 1, 1, checks_per_round=checks, rounds=rounds)
    r = cc.simulate(shots, seed=seed)
    r.update({"d": d, "D": D, "p": p, "seed": seed, "rounds": rounds, "checks": checks, "qudits": cc.c.nq,
              "noise_locations": cc.c.noise_locations()})
    return r


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--shots", type=int, default=2_000_000)
    ap.add_argument("--chunks", type=int, default=10)
    ap.add_argument("--dims", default="2,4,8")
    ap.add_argument("--Ds", default="3")
    ap.add_argument("--ps", default="0.0005,0.001,0.002")
    ap.add_argument("--rounds", type=int, default=2)
    ap.add_argument("--checks", type=int, default=2)
    ap.add_argument("--out", default="runs/cultivation/sweep.json")
    a = ap.parse_args(argv)
    from joblib import Parallel, delayed
    jobs = [(int(d), int(D), float(p), a.shots // a.chunks, 1000 * k + 7, a.rounds, a.checks)
            for d in a.dims.split(",") for D in a.Ds.split(",") for p in a.ps.split(",") for k in range(a.chunks)]
    t0 = time.time()
    res = Parallel(n_jobs=a.jobs, backend="loky")(delayed(one)(*j) for j in jobs)
    agg = {}
    for r in res:
        key = (r["d"], r["D"], r["p"])
        g = agg.setdefault(key, {"d": r["d"], "D": r["D"], "p": r["p"], "shots": 0, "accepted": 0, "logical_failures": 0,
                                 "qudits": r["qudits"], "noise_locations": r["noise_locations"], "rounds": r["rounds"], "checks": r["checks"]})
        for k in ("shots", "accepted", "logical_failures"):
            g[k] += r[k]
    out = []
    for g in agg.values():
        g["accept_rate"] = g["accepted"] / g["shots"]
        g["p_L"] = g["logical_failures"] / max(g["accepted"], 1)
        out.append(g)
        print(f"d={g['d']} D={g['D']} p={g['p']:.4g}: accept {g['accept_rate']:.4f}  p_L {g['p_L']:.3e} ({g['logical_failures']} fails / {g['accepted']})", flush=True)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as f:
        json.dump({"seconds": round(time.time() - t0, 1), "results": out}, f, indent=1)


if __name__ == "__main__":
    main()
