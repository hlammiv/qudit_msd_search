"""p = 3, m = 7, r = 4, linear charge x7 + c (X = two parallel hyperplanes, 1458 points, rank 246):
random output sets at k = 196 (model window: n/k = 6.44, needs d >= 7). Certify weight <= 4
(serial MITM) and weight 5 (parallel 2+3); report d or '>= 6'."""
import sys, random, time, numpy as np
from qmsd.reedmuller import points
from qmsd.weighted import weighted_support, class_points
from qmsd.structured_ad import _rank
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import weight_d_exists_fast, weight5_exists_parallel
p, m, r, k = 3, 7, 4, int(sys.argv[1]); nseeds = int(sys.argv[2])
P = np.asarray(points(m, p), dtype=np.int64); g = (P[:, 6] + 1) % p
GX, X, mu = weighted_support(p, m, g, r); cand = class_points(mu, p, 1)
print(f"|X|={len(X)} rank={_rank(GX, p)} class={len(cand)}", flush=True)
rng = random.Random(3)
for s in range(nseeds):
    cols = sorted(rng.sample(cand, k))
    while _rank(GX[:, cols].T.copy(), p) < k: cols = sorted(rng.sample(cand, k))
    G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p; t0 = time.time(); d = None
    if (G0 == 0).all(axis=0).any(): d = 1
    for dd in (2, 3, 4):
        if d is None and weight_d_exists_fast(G0, p, dd): d = dd
    if d is None: d = 5 if weight5_exists_parallel(G0, p) else 6
    print(f"seed {s}: k={k} n={G0.shape[1]} checks={G0.shape[0]} d={'>=6' if d == 6 else d} gamma(d=7)={np.log(G0.shape[1]/k)/np.log(7):.3f} [{time.time()-t0:.0f}s]", flush=True)
