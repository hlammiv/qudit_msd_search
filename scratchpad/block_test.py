"""Block weight <= 5 check on the sweep's draws (same RNG as plane_sweep.py), reporting the
lightest block word and whether it is planar. usage: block_test.py K SEED0 NSEEDS"""
import sys, random, time, itertools, numpy as np
from qmsd.reedmuller import points
from qmsd.weighted import weighted_support, quadric, class_points
from qmsd.structured_ad import _rank, mitm_streamed
from qmsd.triorthogonal import shorten_matrix
k, seed0, nseeds = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
p, m, r = 7, 4, 7
P = np.asarray(points(m, p), dtype=np.int64)
g = quadric(p, m, (1, 1, 1, 1)); GX, X, mu = weighted_support(p, m, g, r); cand = class_points(mu, p, 1)
for seed in range(seed0, seed0 + nseeds):
    rng = random.Random(1000 * k + seed)
    cols = sorted(rng.sample(cand, k))
    while _rank(GX[:, cols].T.copy(), p) < k: cols = sorted(rng.sample(cand, k))
    G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p
    keep = [l for l in range(GX.shape[1]) if l not in set(cols)]
    cls = [j for j, l in enumerate(keep) if l in set(cand)]
    Gc = G0[:, cls]; Gc = Gc[[t for t in range(Gc.shape[0]) if Gc[t].any()]]
    t0 = time.time(); found = None
    for d in (2, 3, 4, 5):
        ws = list(itertools.islice(mitm_streamed(Gc, d, p), 1))
        if ws:
            supp, vals = ws[0]; pts = [int(X[keep[cls[c]]]) for c in supp]; Q = P[pts]
            found = (d, [tuple(map(int, q)) for q in Q], int(_rank((Q[1:] - Q[0]) % p, p)))
            break
    if found:
        d, qs, span = found
        print(f"k={k} seed={seed}: block weight {d} word, affine span {span} ({'planar' if span <= 2 else 'NON-planar'}) {qs} [{time.time()-t0:.0f}s]", flush=True)
    else:
        print(f"k={k} seed={seed}: no block word of weight <= 5 [{time.time()-t0:.0f}s]", flush=True)
