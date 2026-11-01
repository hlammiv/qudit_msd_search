"""Part 2 as specified: p = 11, m = 4, r = 13 candidates at k = 2294 (26 checks) and k = 2296
(24 checks); build G0, run the Stern finder for weight <= 6 with a fixed iteration budget, record
the model's per-candidate coverage, and try the repo's exact weight-4 scan on two of them.
usage: law_part2.py CANDDIR OUTJSON [iters] [--exact4 N]"""
import sys, os, json, glob, time, numpy as np
from qmsd.reedmuller import rm_generator
from qmsd.triorthogonal import shorten_matrix
from qmsd.structured_ad import _rank
from isd_fp import find_low_weight, expected_words

canddir, outf = sys.argv[1], sys.argv[2]
iters = int(sys.argv[3]) if len(sys.argv) > 3 and not sys.argv[3].startswith("--") else 10
n_exact4 = int(sys.argv[sys.argv.index("--exact4") + 1]) if "--exact4" in sys.argv else 2
state = json.load(open(outf)) if os.path.exists(outf) else {"iters": iters, "candidates": {}}
G_cache = {}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


files = sorted(glob.glob(os.path.join(canddir, "p11_r13_k*_seed*.json")))
exact_done = 0
for f in files:
    name = os.path.basename(f)[:-5]
    c = json.load(open(f)); p, m, r = c["p"], c["m"], c["r"]
    rec = state["candidates"].get(name, {})
    npz = os.path.join(os.path.dirname(outf), name + ".npz")
    if not os.path.exists(npz):
        t = time.time()
        if (r, m, p) not in G_cache:
            G_cache[(r, m, p)] = np.asarray(rm_generator(r, m, p), dtype=np.int64) % p
        G = G_cache[(r, m, p)]
        cols0 = [x - 1 for x in c["puncture_columns_1indexed"]]
        G0 = np.ascontiguousarray(np.asarray(shorten_matrix(G, cols0, p), dtype=np.int64) % p)
        G0 = G0[[t_ for t_ in range(G0.shape[0]) if G0[t_].any()]]
        rk = _rank(G0.copy(), p)
        np.savez_compressed(npz, H=G0, p=p, rank=rk, source=os.path.abspath(f), k=c["k"])
        log(f"{name}: G0 {G0.shape} rank {rk} ({time.time()-t:.0f}s)")
    d = np.load(npz); H = d["H"]; rk = int(d["rank"]); n = H.shape[1]
    rec.update({"k": c["k"], "n": int(n), "checks": rk, "max_outputs_per_plane": c["max_outputs_per_plane"],
                "E_law": {str(w): expected_words(n, rk, p, w) for w in (4, 5, 6, 7)},
                "P_no_w6_law": float(np.exp(-expected_words(n, rk, p, 6)))})
    if "stern6" not in rec:
        res = find_low_weight(H, p, 6, iters=iters, seed=1, method="stern", checks=rk)
        rec["stern6"] = {k: res.get(k) for k in ("found", "support", "weight", "iterations", "seconds", "params", "p_success", "p_miss_one_word", "p_all_missed_model")}
        log(f"{name}: checks {rk} E6 {rec['E_law']['6']:.3g} P(no w6) {rec['P_no_w6_law']:.2f} | stern w<=6 {iters} it: found {res['found']} ({res['seconds']}s) P/iter/word {res['p_success']['6']:.3g}")
    if exact_done < n_exact4 and "exact4" not in rec:
        t = time.time()
        try:
            from qmsd.mindist_nb import weight_d_exists_fast
            ex = {}
            for dd in (2, 3, 4):
                ex[str(dd)] = bool(weight_d_exists_fast(H, p, dd))
                if ex[str(dd)]: break
            rec["exact4"] = {"result": ex, "seconds": round(time.time() - t)}
            log(f"{name}: exact weight<=4: {ex} ({time.time()-t:.0f}s)")
        except MemoryError as e:
            rec["exact4"] = {"error": "MemoryError under ulimit -v 4 GB: " + str(e)[:200], "seconds": round(time.time() - t)}
            log(f"{name}: exact weight-4 scan FAILED (MemoryError) after {time.time()-t:.0f}s")
        except Exception as e:  # noqa
            rec["exact4"] = {"error": type(e).__name__ + ": " + str(e)[:200], "seconds": round(time.time() - t)}
            log(f"{name}: exact weight-4 scan FAILED ({type(e).__name__}: {str(e)[:100]})")
        exact_done += 1
    state["candidates"][name] = rec
    json.dump(state, open(outf, "w"), indent=1, default=float)

# summary table
rows = {}
for name, rec in state["candidates"].items():
    if "stern6" not in rec: continue
    g = rows.setdefault(rec["k"], {"n": 0, "found": 0, "P_no_w6": [], "E6": [], "cov": []})
    g["n"] += 1; g["found"] += int(rec["stern6"]["found"]); g["P_no_w6"].append(rec["P_no_w6_law"]); g["E6"].append(rec["E_law"]["6"])
    g["cov"].append(1 - rec["stern6"]["p_all_missed_model"]["6"])
state["summary"] = {str(k): {"candidates": g["n"], "stern_found_w6": g["found"], "model_P_no_w6": float(np.mean(g["P_no_w6"])),
                             "model_E6": float(np.mean(g["E6"])), "model_P_stern_finds_one": float(np.mean(g["cov"]))} for k, g in rows.items()}
json.dump(state, open(outf, "w"), indent=1, default=float)
print(json.dumps(state["summary"], indent=1))
