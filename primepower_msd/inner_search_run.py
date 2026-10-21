"""Driver: annealing over signed flat collections; records every (n, k, D>=5, normal) code seen per k.

    PYTHONPATH=. python -m primepower_msd.inner_search_run --m 6 --seed 0 --steps 20000 --out runs/meas_distill/inner/m6_s0.json
"""
import argparse
import json
import numpy as np

from .inner_search import anneal, evaluate, signed_rm_flats, FlatPool, compatible


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, default=6)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--steps", type=int, default=20000)
    ap.add_argument("--kmin", type=int, default=16)
    ap.add_argument("--kmax", type=int, default=27)
    ap.add_argument("--restarts", type=int, default=5)
    ap.add_argument("--Dmin", type=int, default=5)
    ap.add_argument("--thin", type=int, nargs=2, default=None, help="start from a random subset of this many RM(2,m) flats")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    rng = np.random.default_rng(a.seed)
    found = {}

    def record(step, e):
        if e and e.get("D", 0) >= a.Dmin and e.get("normal") and a.kmin <= e["k"] <= a.kmax:
            k = e["k"]
            if k not in found or e["n"] < found[k]["n"]:
                found[k] = {kk: e[kk] for kk in ("n", "k", "r", "D", "support", "rows", "ell")}
                found[k]["T_per_out"] = 1 + 2 * ((a.Dmin - 1) // 2 if a.Dmin >= 5 else 2) * e["n"] / k
                print("k", k, "n", e["n"], "D", e["D"], "T/out", round(found[k]["T_per_out"], 3), flush=True)

    for rs in range(a.restarts):
        if a.m == 6 or rs % 2 == 0:
            start, ell = signed_rm_flats(2, a.m)
        else:
            start, ell = signed_rm_flats(2, a.m)
        if a.thin:
            idx = rng.permutation(len(start))[: int(rng.integers(a.thin[0], a.thin[1] + 1))]
            start = [start[i] for i in idx]
        elif a.m >= 7 and a.Dmin <= 5:          # thin out the big RM(2,7) start towards k ~ 20
            idx = rng.permutation(len(start))[: int(rng.integers(12, 22))]
            start = [start[i] for i in idx]
        best = anneal(a.m, start, ell, steps=a.steps, k_range=(a.kmin, a.kmax), seed=int(rng.integers(1 << 30)),
                      log=record, D_min=a.Dmin, s_out=(a.Dmin - 1) // 2 if a.Dmin >= 5 else 2)
        json.dump({"m": a.m, "seed": a.seed, "found": found}, open(a.out, "w"))


if __name__ == "__main__":
    main()
