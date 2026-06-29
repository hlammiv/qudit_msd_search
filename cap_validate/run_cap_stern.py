"""Time-boxed DEEP local Stern/ISD pass on the cap code G0^perp (Method 1).

Sweeps p in {1,2,3,4}; p=3,4 use a restricted (randomised-per-iteration) info_width
to keep the birthday sublists tractable.  Multithreaded across all cores (numba nogil).
Writes a JSON result with the global lowest weight + witness support and per-p history.

Run:  python cap_validate/run_cap_stern.py --minutes 22 --threads 20
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.getcwd())
from cap_validate.stern_isd import build_cap_g0, stern_search, trivial_low_weight  # noqa


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--minutes", type=float, default=22.0)
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 20)
    ap.add_argument("--budget", type=int, default=13)
    ap.add_argument("--out", default="cap_validate/cap_stern_result.json")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    G0 = build_cap_g0()
    r, n = G0.shape
    print(f"G0 {G0.shape}  (r={r}, n={n}, k={n-r})", flush=True)
    triv = trivial_low_weight(G0)
    print(f"trivial: zero_cols={triv['num_zero_columns']} prop_pair="
          f"{triv['proportional_pair']} trivial_min_weight={triv['trivial_min_weight']}",
          flush=True)

    # (p, l, info_width, iters-per-round) -- rotate rounds until the clock runs out.
    # p<=2 enumerate the full k=1913 info set; p=3,4 use a random width window each iter.
    configs = [
        dict(p=1, l=20, info_width=0,   iters=400),
        dict(p=2, l=20, info_width=0,   iters=400),
        dict(p=3, l=18, info_width=260, iters=600),
        dict(p=4, l=16, info_width=80,  iters=800),
    ]

    deadline = time.time() + args.minutes * 60.0
    global_best = args.budget + 1
    global_wit = None
    history = []
    rounds = 0
    per_p_iters = {c["p"]: 0 for c in configs}
    per_p_best = {c["p"]: None for c in configs}
    seed = args.seed

    while time.time() < deadline:
        for c in configs:
            if time.time() >= deadline:
                break
            seed += c["iters"]
            res = stern_search(
                G0, p=c["p"], l=c["l"], iterations=c["iters"],
                weight_budget=args.budget, threads=args.threads, seed=seed,
                target=10, info_width=c["info_width"], verbose=False,
            )
            per_p_iters[c["p"]] += res.iterations
            if res.best_weight != -1:
                if per_p_best[c["p"]] is None or res.best_weight < per_p_best[c["p"]]:
                    per_p_best[c["p"]] = res.best_weight
                if res.best_weight < global_best:
                    global_best = res.best_weight
                    global_wit = res.witness.copy()
                    supp = np.nonzero(global_wit)[0].tolist()
                    history.append(dict(p=c["p"], weight=int(global_best),
                                        round=rounds, support=supp))
                    print(f"[round {rounds} p={c['p']}] NEW BEST weight={global_best} "
                          f"supp={supp}", flush=True)
                    # verify on the fly
                    syn = (G0.astype(np.int64) @ (global_wit.astype(np.int64) % 3)) % 3
                    assert np.count_nonzero(syn) == 0, "witness not a codeword!"
                    if global_best < 10:
                        print("  *** weight < 10 -> d < 10 REFUTED ***", flush=True)
        rounds += 1
        elapsed = time.time() - (deadline - args.minutes * 60.0)
        tot = sum(per_p_iters.values())
        print(f"  ... round {rounds} done, total_iters={tot}, "
              f"best={global_best if global_best <= args.budget else None}, "
              f"{elapsed:.0f}s", flush=True)
        if global_best < 10:
            break

    bw = int(global_best) if global_best <= args.budget else -1
    supp = np.nonzero(global_wit)[0].tolist() if global_wit is not None else None
    total_iters = sum(per_p_iters.values())
    result = dict(
        method="stern_isd",
        code="[[1968,219]]_3 (m=7 cap qutrit)",
        n=n, r=r, k=n - r,
        threads=args.threads,
        minutes=args.minutes,
        weight_budget=args.budget,
        lowest_weight_found=bw,
        below_10=(bw != -1 and bw < 10),
        witness_support=supp,
        total_iterations=total_iters,
        rounds=rounds,
        per_p_iterations=per_p_iters,
        per_p_best=per_p_best,
        history=history,
        elapsed_sec=time.time() - (deadline - args.minutes * 60.0),
        verdict=("refuted-d-below-10" if (bw != -1 and bw < 10)
                 else "no-counterexample-yet"),
    )
    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps({k: v for k, v in result.items() if k != "history"}, indent=2,
                     default=int), flush=True)
    print(f"WROTE {args.out}", flush=True)


if __name__ == "__main__":
    main()
