"""r = 6, g = product of 5 general hyperplanes (|X| = 1110, rank 205): random cube-class output
sets at k = 180, 190; certify weight <= 4 (serial) and weight 5 (parallel full 2+3 scan)."""
import random, time, json, numpy as np
from qmsd.weighted import weighted_support, product_of_linear, class_points
from qmsd.structured_ad import _rank
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import weight_d_exists_fast, weight5_exists_parallel
p, m, r = 7, 4, 6
g = product_of_linear(p, m, [((1, 0, 0, 0), 0), ((0, 1, 0, 0), 0), ((0, 0, 1, 0), 0), ((1, 1, 1, 1), 2), ((1, 2, 4, 1), 3)])
GX, X, mu = weighted_support(p, m, g, r); cand = class_points(mu, p, 1)
print(f"|X|={len(X)} rank={_rank(GX, p)} class pts={len(cand)} rank on class={_rank(GX[:, cand], p)}", flush=True)
rng = random.Random(11); out = []
for k in (180, 190):
    for s in range(2):
        cols = sorted(rng.sample(cand, k))
        while _rank(GX[:, cols].T.copy(), p) < k: cols = sorted(rng.sample(cand, k))
        G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p; n = G0.shape[1]
        t0 = time.time(); d = None
        for dd in (2, 3, 4):
            if weight_d_exists_fast(G0, p, dd): d = dd; break
        if d is None:
            d = 5 if weight5_exists_parallel(G0, p) else 6
        print(f"k={k} s={s}: n={n} dimG0={G0.shape[0]} d={'>=6' if d == 6 else d} gamma={np.log(n/k)/np.log(d):.3f} [{time.time()-t0:.0f}s]", flush=True)
        out.append({"k": k, "sample": s, "n": n, "d": d, "cols_in_X": cols}); json.dump(out, open("runs/weighted/p7_m4_r6_hyperplanes_probe.json", "w"), indent=1)
