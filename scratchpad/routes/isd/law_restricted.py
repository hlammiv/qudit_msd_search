"""Exact law test at p = 11, m = 4 in a regime where the exhaustive scan is feasible: the
RM(13, 4) kernel code of a k = 2294 candidate is shortened further (extra outputs) down to
``checks`` rows and restricted to random column subsets of size s; on each subset the weight-6
words are counted EXACTLY by the 3+3 projective MITM (after excluding weight <= 5).  The law
predicts Poisson(E_6) counts with E_6 = C(s, 6) (p-1)^5 / p^checks (checks = rank of the
restriction).  Same algebraic family as the m = 4 candidates (evaluation vectors of kernel
polynomials at random points of F_11^4); structured (plane) words have weight >= 55 - occupancy
so cannot contaminate weight 6.
usage: law_restricted.py CAND.npz OUTJSON checks s nsubsets [seed]"""
import sys, os, json, time, numpy as np
from qmsd.triorthogonal import shorten_matrix
from qmsd.structured_ad import _rank
from isd_fp import mitm_exhaustive, expected_words, list_size

npz, outf = sys.argv[1], sys.argv[2]
checks, s, nsub = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
seed = int(sys.argv[6]) if len(sys.argv) > 6 else 0
d = np.load(npz); H = d["H"]; p = int(d["p"]); r0 = int(d["rank"]); n = H.shape[1]
rng = np.random.default_rng(seed)
# extra shortening: kill r0 - checks random columns (they become outputs)
extra = np.sort(rng.choice(n, r0 - checks, replace=False))
Hs = np.asarray(shorten_matrix(H, [int(x) for x in extra], p), dtype=np.int64) % p
Hs = Hs[[t for t in range(Hs.shape[0]) if Hs[t].any()]]
rest = np.array([c for c in range(n) if c not in set(extra.tolist())])
assert Hs.shape[1] == len(rest)
rk = _rank(Hs.copy(), p)
assert rk == checks, (rk, checks)
state = json.load(open(outf)) if os.path.exists(outf) else {"source": npz, "p": p, "checks": checks, "s": s, "extra_outputs": extra.tolist(), "subsets": []}
E = {w: expected_words(s, checks, p, w) for w in (4, 5, 6, 7)}
print(f"p={p} checks={checks} s={s}: law E4={E[4]:.3g} E5={E[5]:.3g} E6={E[6]:.3g} E7={E[7]:.3g}; 3+3 list {list_size(s, 3, p):.3g}", flush=True)
state["E_law"] = {str(w): v for w, v in E.items()}
done = len(state["subsets"])
sub_rng = np.random.default_rng(seed + 1000)
for i in range(nsub):
    S = np.sort(sub_rng.choice(len(rest), s, replace=False))
    if i < done: continue
    cols = S                                   # positions in Hs (shorten_matrix dropped the extra columns)
    Hsub = np.ascontiguousarray(Hs[:, cols])
    rks = _rank(Hsub.copy(), p)
    t = time.time(); low = None; words = []
    for w in (2, 3, 4, 5):
        ws = mitm_exhaustive(Hsub, p, w, rng=rng, return_all=True)
        if ws: low = w; words = ws; break
    if low is None:
        words = mitm_exhaustive(Hsub, p, 6, rng=rng, return_all=True)
    rec = {"subset_index": i, "rank": int(rks), "lowest_weight": low or (6 if words else None), "n_words": len(words),
           "words": [[int(x) for x in w_[0]] for w_ in words], "seconds": round(time.time() - t, 1)}
    state["subsets"].append(rec)
    json.dump(state, open(outf, "w"), indent=1)
    print(f"  subset {i}: rank {rks} lowest weight {rec['lowest_weight']} count {len(words)} ({rec['seconds']}s)", flush=True)
counts = [r_["n_words"] if r_["lowest_weight"] in (None, 6) else -1 for r_ in state["subsets"]]
c6 = [c for c in counts if c >= 0]
state["summary"] = {"subsets": len(counts), "with_weight_le5": sum(1 for c in counts if c < 0), "weight6_counts": c6,
                    "mean_count": float(np.mean(c6)) if c6 else None, "frac_zero": float(np.mean([c == 0 for c in c6])) if c6 else None,
                    "law_E6": E[6], "law_P_zero": float(np.exp(-E[6]))}
json.dump(state, open(outf, "w"), indent=1)
print(json.dumps(state["summary"]))
