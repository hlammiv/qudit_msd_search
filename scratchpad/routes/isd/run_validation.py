"""Part 1 validation runs of the ISD finder on the three known instances (see NOTES.md).
usage: run_validation.py OUTDIR  [which: p7a p7b p11 ...]"""
import sys, os, json, time, numpy as np
from isd_fp import (find_low_weight, stern_iteration, expected_words, p_success_stern, cost_table,
                    list_size, verify_word)

outdir = sys.argv[1]; which = set(sys.argv[2:]) or {"p7a", "p7b", "p11"}
RUNS = os.path.join(outdir)


def load(name):
    d = np.load(os.path.join(RUNS, name)); return d["H"], int(d["p"]), int(d["rank"])


def save(name, obj):
    json.dump(obj, open(os.path.join(RUNS, name), "w"), indent=1, default=float)
    print("saved", name, flush=True)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ---- (a) p = 7, k = 272 seed 3 class block: a weight-6 word exists at [16, 67, 122, 200, 210, 328]
if "p7a" in which:
    H, p, r = load("p7_k272_seed3_block.npz"); n = H.shape[1]
    known = [16, 67, 122, 200, 210, 328]
    rk = np.linalg.matrix_rank(H[:, known].astype(float))  # just informative; exact check below
    out = {"instance": "p7_k272_seed3_block", "n": n, "checks": r, "known_word_cols": known,
           "E_law": {str(w): expected_words(n, r, p, w) for w in (5, 6, 7)}, "runs": []}
    # exact: is the known support really dependent?
    from qmsd.structured_ad import _null_space
    ns = _null_space(H[:, known], p); out["known_word_nullity"] = int(ns.shape[0])
    log("p7a: n", n, "checks", r, "E6", out["E_law"]["6"], "known support nullity", ns.shape[0])
    for seed in range(5):
        res = find_low_weight(H, p, 6, iters=5000, seed=seed, method="stern", checks=r)
        log(f"p7a stern seed {seed}: found {res['found']} {res['support']} iters {res['iterations']} {res['seconds']}s")
        out["runs"].append({k: res[k] for k in ("found", "support", "coefficients", "weight", "iterations", "seconds", "params", "p_success")})
    # how many DISTINCT weight-6 words does the block have? Stern with --all for a fixed budget
    rng = np.random.default_rng(100); words = {}; t = time.time(); its = 0
    while time.time() - t < 240:
        ws, info = stern_iteration(H, p, 6, rng, 2, 2); its += 1
        for s, c in ws: words.setdefault(tuple(s), 0); words[tuple(s)] += 1
    out["census_240s"] = {"iterations": its, "distinct_weight6_words": sorted([list(k) for k in words]),
                          "hits_per_word": {str(k): v for k, v in words.items()},
                          "p_success_6": p_success_stern(n, r, 6, 2, 2, info["ell"], info["z1"], info["z2"])}
    log(f"p7a census: {its} iters, distinct weight-6 words {len(words)}: {list(words.items())}")
    save("p7a_k272_seed3_stern.json", out)

# ---- (b) p = 7, k = 270 seed 19 class block: no word of weight <= 6, weight-7 words exist (~27 predicted)
if "p7b" in which:
    H, p, r = load("p7_k270_seed19_block.npz"); n = H.shape[1]
    out = {"instance": "p7_k270_seed19_block", "n": n, "checks": r,
           "E_law": {str(w): expected_words(n, r, p, w) for w in (5, 6, 7, 8)}, "runs": []}
    log("p7b: n", n, "checks", r, "E6", out["E_law"]["6"], "E7", out["E_law"]["7"])
    res = find_low_weight(H, p, 6, iters=3000, seed=0, method="stern", checks=r)
    log(f"p7b stern w<=6: found {res['found']} iters {res['iterations']} {res['seconds']}s p_success {res['p_success']} p_miss_one {res.get('p_miss_one_word')}")
    out["w6"] = {k: res.get(k) for k in ("found", "support", "iterations", "seconds", "params", "p_success", "p_miss_one_word", "p_all_missed_model")}
    # weight 7: first hit, then a fixed-budget census of distinct words (Monte-Carlo estimate of N_7)
    res = find_low_weight(H, p, 7, iters=5000, seed=1, method="stern", checks=r)
    log(f"p7b stern w<=7: found {res['found']} {res['support']} iters {res['iterations']} {res['seconds']}s p_success {res['p_success']}")
    out["w7_first"] = {k: res.get(k) for k in ("found", "support", "coefficients", "weight", "iterations", "seconds", "params", "p_success")}
    rng = np.random.default_rng(7); words = {}; t = time.time(); its = 0; hits = 0
    while its < 3000:
        ws, info = stern_iteration(H, p, 7, rng, 2, 2); its += 1
        for s, c in ws:
            words.setdefault(tuple(s), 0); words[tuple(s)] += 1; hits += 1
    P7 = p_success_stern(n, r, 7, 2, 2, info["ell"], info["z1"], info["z2"])
    P6 = p_success_stern(n, r, 6, 2, 2, info["ell"], info["z1"], info["z2"])
    w6_in_census = [list(k) for k in words if len(k) <= 6]
    out["w7_census"] = {"iterations": its, "seconds": round(time.time() - t), "total_hits": hits, "distinct_words": len(words),
                        "weights": sorted(len(k) for k in words), "weight_le6_words": w6_in_census,
                        "p_success_7": P7, "p_success_6": P6,
                        "N7_estimate_hits_over_itsP": hits / (its * P7), "E7_law": out["E_law"]["7"],
                        "words": {str(list(k)): v for k, v in words.items()}}
    log(f"p7b w7 census: {its} iters {hits} hits, {len(words)} distinct; N7 estimate {hits/(its*P7):.1f} vs law {out['E_law']['7']:.1f}; weight<=6 seen: {w6_in_census}")
    save("p7b_k270_seed19_stern.json", out)

# ---- (c) p = 11, G0 30 x 12351 from runs/rm_m4/p11_r13_k2290_seed1.json
if "p11" in which:
    H, p, r = load("p11_r13_k2290_seed1.npz"); n = H.shape[1]
    out = {"instance": "p11_r13_k2290_seed1", "n": n, "checks": r,
           "E_law": {str(w): expected_words(n, r, p, w) for w in (4, 5, 6, 7, 8)}, "cost_table": {}, "runs": {}}
    for w in (6, 7):
        out["cost_table"][str(w)] = [{"variant": a, "list": b, "p_per_iter_per_word": c} for a, b, c in cost_table(n, r, p, w)]
    log("p11: n", n, "checks", r, "E", out["E_law"])
    for w, iters in ((6, 40), (7, 40)):
        res = find_low_weight(H, p, w, iters=iters, seed=0, method="stern", checks=r)
        log(f"p11 stern w<={w}: found {res['found']} iters {res['iterations']} {res['seconds']}s p_success {res['p_success']} p_all_missed_model {res.get('p_all_missed_model')}")
        out["runs"][str(w)] = {k: res.get(k) for k in ("found", "support", "iterations", "seconds", "params", "p_success", "p_miss_one_word", "p_all_missed_model")}
    save("p11_k2290_stern.json", out)
