"""Lenore: scale plain-cap sampling to collect the larger paper qutrit codes that plane_spread
is too restrictive for (plane-spread caps only exist to k~17).
Targets: [[206,37,4]]_3 (k=37, capset reaches d=4 ~1/640) and [[215,28,5]]_3 (k=28, d=5 rare
needle among caps). At p=3,m=5 the MITM is small -> full 28-way parallelism, no memory issue.
Records the puncture columns of any code reaching the paper distance. Incremental JSON save."""
import os, json, time, math, random
import numpy as np
from joblib import Parallel, delayed
from qmsd.reedmuller import r_max, rm_generator
from qmsd.codes import code_from_puncture
from qmsd.sampling import random_cap, points_to_columns, all_points

p, m = 3, 5
r = r_max(m, p); N = p ** m
G = rm_generator(r, m, p)
allpts = all_points(m, p)
TARGETS = [(37, 4, "[[206,37,4]]"), (28, 5, "[[215,28,5]]")]
NJOBS = 28
TIME_BUDGET = 2.0 * 3600
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "capset_scale_results.json")

def one(k, seed, want_cols):
    pts = random_cap(m, p, k, random.Random(seed), allpts)
    if pts is None:
        return (k, -1, None)                      # cap stall
    cols = points_to_columns(pts, p)
    c = code_from_puncture(p, m, cols, r=r, compute_A_d=False, max_distance=6, G=G)
    if not (c.full_rank and c.d):
        return (k, 0, None)
    return (k, int(c.d), list(cols) if c.d >= want_cols else None)

def main():
    state = {"meta": {"p": p, "m": m, "started": time.time(), "draws": 0},
             "by_k": {str(k): {"target": dp, "label": lbl, "built": 0, "stalls": 0, "max_d": 0, "hist": {}}
                      for k, dp, lbl in TARGETS},
             "found": []}
    t0 = time.time(); batch = 0
    while time.time() - t0 < TIME_BUDGET:
        jobs = []
        for k, dp, _ in TARGETS:
            for j in range(NJOBS):
                jobs.append((k, batch * 100003 + j * 7 + k, dp))
        out = Parallel(n_jobs=NJOBS)(delayed(one)(k, s, dp) for (k, s, dp) in jobs)
        for (k, d, cols) in out:
            bk = state["by_k"][str(k)]
            state["meta"]["draws"] += 1
            if d == -1:
                bk["stalls"] += 1; continue
            if d == 0:
                continue
            bk["built"] += 1
            bk["hist"][str(d)] = bk["hist"].get(str(d), 0) + 1
            bk["max_d"] = max(bk["max_d"], d)
            if cols is not None:
                g = math.log((N - k) / k) / math.log(d)
                state["found"].append({"label": bk["label"], "n": N - k, "k": k, "d": d,
                                        "gamma": round(g, 4), "cols": cols})
        batch += 1
        state["meta"]["elapsed_min"] = round((time.time() - t0) / 60, 1)
        with open(OUT, "w") as f:
            json.dump(state, f)
        msg = "  ".join(f"{state['by_k'][str(k)]['label']}: max_d={state['by_k'][str(k)]['max_d']} "
                        f"(built {state['by_k'][str(k)]['built']}, stalls {state['by_k'][str(k)]['stalls']})"
                        for k, _, _ in TARGETS)
        print(f"[{state['meta']['elapsed_min']:6.1f}m] {msg}  found={len(state['found'])}", flush=True)
    print(f"DONE. codes reaching paper distance: {len(state['found'])}", flush=True)

if __name__ == "__main__":
    main()
