"""Lenore: scaled NEAR-cap hunt for [[200,43,3]]_3 (k=43, d=3 rare needle).
The paper's set is a near-cap with 31 collinear triples; strict caps stall at k=43. Near-caps
(triple budget t) build fine but give d=2 -> d=3 is rare (0/120 local). Scale it. Parallel,
incremental save; record the puncture columns of any d=3 code found."""
import os, json, time, math, random
from joblib import Parallel, delayed
from qmsd.reedmuller import r_max, rm_generator
from qmsd.codes import code_from_puncture
from qmsd.sampling import all_points, collinear, points_to_columns

p, m = 3, 5
r = r_max(m, p); N = p ** m
G = rm_generator(r, m, p)
allpts = all_points(m, p)
K = 43
TARGET_D = 3
TVALS = [31, 36, 42]          # triple budgets (paper's set has 31)
NJOBS = 28
TIME_BUDGET = 1.5 * 3600
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nearcap_43_results.json")

def new_coll(chosen, x):
    c = 0
    for i in range(len(chosen)):
        a = chosen[i]
        for j in range(i + 1, len(chosen)):
            if collinear(a, chosen[j], x, p):
                c += 1
    return c

def near_cap(k, t, rng):
    order = list(allpts); rng.shuffle(order); chosen = []; used = 0
    for x in order:
        nc = new_coll(chosen, x)
        if used + nc <= t:
            chosen.append(x); used += nc
            if len(chosen) == k:
                return chosen
    return None

def one(t, seed):
    pts = near_cap(K, t, random.Random(seed))
    if pts is None:
        return (t, -1, None)                       # stall
    cols = points_to_columns(pts, p)
    c = code_from_puncture(p, m, cols, r=r, compute_A_d=False, max_distance=4, G=G)
    if not (c.full_rank and c.d):
        return (t, 0, None)
    return (t, int(c.d), list(cols) if c.d >= TARGET_D else None)

def main():
    state = {"meta": {"p": p, "m": m, "k": K, "started": time.time(), "draws": 0},
             "by_t": {str(t): {"built": 0, "stalls": 0, "hist": {}} for t in TVALS},
             "found": []}
    t0 = time.time(); batch = 0
    while time.time() - t0 < TIME_BUDGET and not state["found"]:   # stop as soon as one d=3 is found
        jobs = [(t, batch * 100003 + j * 13 + t) for t in TVALS for j in range(NJOBS)]
        out = Parallel(n_jobs=NJOBS)(delayed(one)(t, s) for (t, s) in jobs)
        for (t, d, cols) in out:
            bt = state["by_t"][str(t)]
            state["meta"]["draws"] += 1
            if d == -1:
                bt["stalls"] += 1; continue
            if d == 0:
                continue
            bt["built"] += 1
            bt["hist"][str(d)] = bt["hist"].get(str(d), 0) + 1
            if cols is not None:
                g = math.log((N - K) / K) / math.log(d)
                state["found"].append({"label": "[[200,43,3]]", "n": N - K, "k": K, "d": d,
                                       "gamma": round(g, 4), "t": t, "cols": cols})
        batch += 1
        state["meta"]["elapsed_min"] = round((time.time() - t0) / 60, 1)
        with open(OUT, "w") as f:
            json.dump(state, f)
        msg = "  ".join(f"t{t}: built {state['by_t'][str(t)]['built']} hist {state['by_t'][str(t)]['hist']} "
                        f"stalls {state['by_t'][str(t)]['stalls']}" for t in TVALS)
        print(f"[{state['meta']['elapsed_min']:5.1f}m] {msg}  FOUND={len(state['found'])}", flush=True)
    if state["found"]:
        f = state["found"][0]
        print(f"*** REPRODUCED [[200,43,3]] d=3 (t={f['t']}, gamma={f['gamma']}) -- 7/7 qutrit codes ***", flush=True)
    else:
        print(f"DONE (budget) -- no d=3 found in {state['meta']['draws']} draws", flush=True)

if __name__ == "__main__":
    main()
