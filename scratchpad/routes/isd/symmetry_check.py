"""Is the weight-6 excess of the k = 272 shift-orbit block (N_6 ~ 115 vs law 3.7) a property of
shift-symmetric output sets or of k = 272 in general?  Monte-Carlo census (Stern 2+2, fixed time
budget) of distinct weight-6 class-block words on: the shift-orbit draws seed 18 and 40 of
runs/weighted/orbit_sweep/orbit_shift_k272.jsonl, and fresh RANDOM full-rank k = 272 class-point
sets.  N_6 is estimated as hits / (iterations * P_success) and by the distinct count.
usage: symmetry_check.py OUTJSON [budget_seconds]"""
import sys, json, time, numpy as np
from qmsd.weighted import weighted_support, quadric, class_points
from qmsd.structured_ad import _rank
from qmsd.triorthogonal import shorten_matrix
from isd_fp import stern_iteration, p_success_stern, expected_words

outf = sys.argv[1]; budget = float(sys.argv[2]) if len(sys.argv) > 2 else 150.0
p, m, r = 7, 4, 7
g = quadric(p, m, (1, 1, 1, 1)); GX, X, mu = weighted_support(p, m, g, r)
cand = class_points(mu, p, 1); cset = set(cand)


def block(cols):
    G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p
    keep = [l for l in range(GX.shape[1]) if l not in set(cols)]
    C_cols = [j for j in range(G0.shape[1]) if keep[j] in cset]
    Gc = np.ascontiguousarray(G0[:, C_cols]); Gc = Gc[[t for t in range(Gc.shape[0]) if Gc[t].any()]]
    return Gc, int(_rank(Gc.copy(), p))


def census(H, rk, label):
    n = H.shape[1]; rng = np.random.default_rng(1); words = {}; t = time.time(); its = 0; hits = 0; info = None
    while time.time() - t < budget:
        ws, info = stern_iteration(H, p, 6, rng, 2, 2); its += 1
        for s, c in ws:
            words[tuple(s)] = words.get(tuple(s), 0) + 1; hits += 1
    P = p_success_stern(n, rk, 6, 2, 2, info["ell"], info["z1"], info["z2"])
    rec = {"label": label, "n": n, "checks": rk, "iterations": its, "hits": hits, "distinct": len(words),
           "N6_estimate": hits / (its * P), "p_success_6": P, "E6_law": expected_words(n, rk, p, 6),
           "min_weight_seen": min((len(k) for k in words), default=None), "words": [list(k) for k in words]}
    print(f"{label}: n {n} checks {rk} | {its} iters, {hits} hits, {len(words)} distinct, N6 est {rec['N6_estimate']:.1f}, law {rec['E6_law']:.2f}", flush=True)
    return rec


out = []
recs = {}
for l in open("runs/weighted/orbit_sweep/orbit_shift_k272.jsonl"):
    if l.strip():
        rr = json.loads(l)
        if "cols_in_X" in rr: recs[rr["seed"]] = rr
for seed in (18, 40, 7, 25):
    if seed not in recs: continue
    Gc, rk = block(recs[seed]["cols_in_X"])
    out.append(census(Gc, rk, f"shift-orbit k=272 seed {seed} (d_plane {recs[seed]['d_plane']}, block_d_le5 {recs[seed]['block_d_le5']})"))
    json.dump(out, open(outf, "w"), indent=1)
rng = np.random.default_rng(2024)
made = 0
while made < 3:
    cols = sorted(int(x) for x in rng.choice(cand, 272, replace=False))
    Gc, rk = block(cols)
    if rk != 291 - 272: continue
    made += 1
    out.append(census(Gc, rk, f"random k=272 class set #{made}"))
    json.dump(out, open(outf, "w"), indent=1)
