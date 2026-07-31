"""Lenore hunt: is there a gamma<1 code in the p=7 m=4 mid-k window (high-k end)?

Distance-MAXIMISING swap climb over puncture sets S (caps of size >=300 do NOT exist in
AG(4,7), so no cap constraint -- just full-rank + distance-increasing swaps), with EXACT
certification by the memory-balanced MITM (n_jobs internal parallelism; handles a d>=7 hit
without OOM, unlike the unbalanced engine's 14GB table). Runs on lenore (125GB, 32 cores).

gamma<1 needs d > n/k: k=304..312 -> needed d = 7 (or 8 at k<=300). We hunt d>=needed.
Any hit is a NEW gamma<1 p=7 search code. Otherwise the max achievable d is closure evidence.
Incremental JSON save after every trial -> survives ssh drops.
"""
import os, json, time, math, random
import numpy as np
from qmsd.reedmuller import r_max, rm_generator
from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.mindist_balanced import min_dependent_columns_balanced

p, m = 7, 4
r = r_max(m, p)
N = p ** m
KS = [312, 311, 310, 309, 308, 307, 306, 305, 304]   # highest-k (cheapest, best gamma<1 chance) first
NJOBS = 28            # internal parallelism of the balanced certifier
DMAX = 6              # d_max=6: a <=6 dependency => d<=6 (gamma>=1); none found => d>=7 (gamma<1!)
CLIMB_STEPS = 25
SWAP_TRIES = 4
TIME_BUDGET = 5.0 * 3600
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "p7_hunt_results.json")

G = rm_generator(r, m, p)

def needed_d(k):
    return math.floor((N - k) / k) + 1

def certify(cols):
    """exact d via balanced MITM; returns int d (2..6) or 7 meaning d>=7 (gamma<1 candidate), or None if not full rank."""
    built = build_triorthogonal_code(p, m, r, cols, G=G)
    if not built["full_rank"]:
        return None
    G0 = np.asarray(built["X_stab"], dtype=int) % p
    d = min_dependent_columns_balanced(G0, p, d_max=DMAX, n_jobs=NJOBS, ram_fraction=0.55)
    if d is None:
        return 7
    return int(d)

def climb(k, seed, deadline, state):
    rng = random.Random(seed)
    cur = tuple(sorted(rng.sample(range(1, N + 1), k)))
    cur_d = certify(cur)
    while cur_d is None:                       # reseed until full rank
        cur = tuple(sorted(rng.sample(range(1, N + 1), k))); cur_d = certify(cur)
    _record(state, k, cur_d, cur)
    for _ in range(CLIMB_STEPS):
        if time.time() > deadline: break
        cur_set = set(cur)
        outside = [c for c in range(1, N + 1) if c not in cur_set]
        rng.shuffle(outside)
        drop = rng.randrange(len(cur))
        base = list(cur[:drop]) + list(cur[drop + 1:])
        moved = False
        for nc in outside[:SWAP_TRIES]:
            cand = tuple(sorted(base + [nc]))
            d = certify(cand)
            if d is not None and d >= cur_d:
                cur, cur_d, moved = cand, d, True
                _record(state, k, cur_d, cur)
                break
        if not moved:
            break
    return cur_d

def _record(state, k, d, cols):
    bk = state["by_k"][str(k)]
    bk["n"] += 1
    bk["hist"][str(d)] = bk["hist"].get(str(d), 0) + 1
    if d > bk["max_d"]:
        bk["max_d"] = d
    if d >= bk["needed_d"]:
        state["candidates"].append({"k": k, "d": d, "n": N - k,
                                    "gamma": math.log((N - k) / k) / math.log(d), "cols": list(cols)})
    state["meta"]["trials"] += 1
    with open(OUT, "w") as f:
        json.dump(state, f)

def main():
    state = {"meta": {"p": p, "m": m, "started": time.time(), "trials": 0}, "by_k": {}, "candidates": []}
    for k in KS:
        state["by_k"][str(k)] = {"needed_d": needed_d(k), "max_d": 0, "n": 0, "hist": {}}
    t0 = time.time(); deadline = t0 + TIME_BUDGET; s = 0
    while time.time() < deadline:
        for k in KS:
            if time.time() >= deadline: break
            best = climb(k, 7000 + s, deadline, state)
            el = (time.time() - t0) / 60
            print(f"[{el:6.1f}m] k={k} needed>={needed_d(k)} climb best d={best}  "
                  f"(k max_d={state['by_k'][str(k)]['max_d']}, n={state['by_k'][str(k)]['n']})  "
                  f"candidates={len(state['candidates'])}", flush=True)
            s += 1
    print(f"DONE elapsed {(time.time()-t0)/3600:.2f}h trials={state['meta']['trials']} "
          f"gamma<1 candidates={len(state['candidates'])}", flush=True)
    for k in KS:
        bk = state["by_k"][str(k)]
        print(f"  k={k}: max_d={bk['max_d']} (needed {bk['needed_d']}) over n={bk['n']}  hist={bk['hist']}", flush=True)

if __name__ == "__main__":
    main()
