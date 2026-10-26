"""Weighted p=7 m=3 candidates: find a d=5 cap, score by the per-output surrogate (one-pass MITM)
and by marginal_saving vs the p=7 standard library (cyclic + RS/RM + doubled [[20,1,4]]_7)."""
import sys, time, json, math
import numpy as np
sys.path.insert(0, ".")
from qmsd.reedmuller import points, r_max, rm_generator
from qmsd.weighted import weighted_support, weighted_code, quadric, class_points
from qmsd.sampling import random_flat_spread
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist import min_dependent_columns
from qmsd.performance import canonical_logical_rows
from qmsd.objective import surrogate_from_matrices, standard_library, marginal_saving
from qmsd.field import GFp
p, m = 7, 3
GF = GFp(p)
Pa = np.asarray(points(m, p), dtype=np.int64)
P = [tuple(int(v) for v in x) for x in Pa]
index = {x: i for i, x in enumerate(P)}
cases = [
    ("ellipsoid x1^2+x2^2+x3^2+1 [2*cube]", quadric(p, m, (1, 1, 1), const=1), 2, 24),
    ("ellipsoid x1^2+x2^2+x3^2+1 [2*cube]", quadric(p, m, (1, 1, 1), const=1), 2, 21),
    ("parabolic x1^2+x2^2+x3 [cube]", quadric(p, m, (1, 1, 0), (0, 0, 1)), 1, 24),
    ("x1(x1-1) [cube]", quadric(p, m, (1, 0, 0), const=5), 1, 10),
    ("unit", np.ones(len(P), dtype=np.int64), 1, 24),
]
base = standard_library(p)
print("standard library:", [c.name for c in base], flush=True)
rng = np.random.default_rng(13)
results = []
for name, g, unit, k in cases:
    GX, X, mu = weighted_support(p, m, g)
    pos = {int(c): i for i, c in enumerate(X)}
    cand = [P[X[i]] for i in class_points(mu, p, unit)]
    found = None
    for trial in range(60):
        pts = random_flat_spread(m, p, k, rng, allpts=cand)
        if pts is None: continue
        cols = [pos[index[x]] for x in pts]
        if np.linalg.matrix_rank(GF(GX[:, cols])) < k: continue
        G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p
        try:
            d = min_dependent_columns(G0, p, d_max=5)
        except ValueError:
            d = 6
        if d >= 5:
            found = (cols, pts, G0, d); break
    if found is None:
        print(f"{name} k={k}: no d=5 cap in 60 trials", flush=True); continue
    cols, pts, G0, d = found
    H1 = np.asarray(canonical_logical_rows(GX, cols, p), dtype=np.int64) % p
    t0 = time.time()
    cur = surrogate_from_matrices(G0, H1, p, d=d, name=f"{name} [[{G0.shape[1]},{k},{d}]]_7")
    t1 = time.time() - t0
    print(f"\n{cur.name}: gamma={cur.gamma:.3f} L_d={cur.L_d:.1f} c_eps={cur.c_eps:.4g} threshold={cur.threshold_eps:.4f} pacc(1e-2)={cur.pacc(1e-2):.3f} [mitm {t1:.0f}s]", flush=True)
    res = marginal_saving(cur, base)
    print(f"  marginal saving best = {res['best']:.3f}", flush=True)
    for r in res["rows"]:
        print(f"    eps0={r['eps0']:g} target={r['target']:g}: base {r['base']:.4g} with {r['with']:.4g} saving {r['saving']:.3f} uses={r['uses']}", flush=True)
    results.append({"charge": name, "k": k, "n": int(G0.shape[1]), "d": int(d), "gamma": cur.gamma, "L_d": cur.L_d,
                    "c_eps": cur.c_eps, "threshold_eps": cur.threshold_eps, "saving": res,
                    "puncture_points": [list(x) for x in pts], "puncture_cols_in_X": [int(c) for c in cols],
                    "X_global_1indexed": [int(x) + 1 for x in X], "mu": [int(u) for u in mu]})
    json.dump(results, open("weighted_score_results.json", "w"), indent=1)
