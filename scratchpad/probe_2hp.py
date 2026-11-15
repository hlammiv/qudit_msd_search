"""Single-species weighted support at p = 7, m = 4, r = 6: g = x1(x1-1)(x1-2)(x1-3)(x1-4) is +1 on x1 = 5
and -1 on x1 = 6 (U and U^dagger, one Clifford species), |X| = 686, rank 140.  Random output sets at
several k; d of ker(G0) by weight <= 4 serial MITM and the weight-5 parallel 2+3 scan.
usage: PYTHONPATH=. python scratchpad/probe_2hp.py K1,K2,.. DRAWS SEED OUT.json"""
import sys, random, time, json, itertools, numpy as np
from qmsd.reedmuller import points
from qmsd.weighted import weighted_support, product_of_linear
from qmsd.structured_ad import _rank
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import weight_d_exists_fast, weight5_exists_parallel

ks = [int(v) for v in sys.argv[1].split(",")]; draws, seed, out_path = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
p, m, r = 7, 4, 6
e0 = (1, 0, 0, 0)
g = product_of_linear(p, m, [(e0, c) for c in range(5)])
GX, X, mu = weighted_support(p, m, g, r)
P = np.asarray(points(m, p), dtype=np.int64)[X]
assert set(int(v) for v in mu) == {1, 6}, set(mu)
# weighted triorthogonality <=> sum_X mu * x^a = 0 for every monomial of degree <= 3r
for a in itertools.product(range(p), repeat=m):
    if sum(a) <= 3 * r:
        s = int(np.sum(mu * np.prod([pow_col for pow_col in (np.array([pow(int(v), ai, p) for v in P[:, i]]) for i, ai in enumerate(a))], axis=0)) % p)
        assert s == 0, a
print(f"|X|={len(X)} rank={_rank(GX, p)} charges +1: {(mu == 1).sum()}, -1: {(mu == 6).sum()}; triorthogonal OK", flush=True)
rng = random.Random(seed); out = []
cand = list(range(len(X)))
for k in ks:
    for s in range(draws):
        cols = sorted(rng.sample(cand, k))
        while _rank(GX[:, cols].T.copy(), p) < k: cols = sorted(rng.sample(cand, k))
        G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p; n = G0.shape[1]
        t0 = time.time(); d = None
        for dd in (2, 3, 4):
            if weight_d_exists_fast(G0, p, dd): d = dd; break
        if d is None:
            d = 5 if weight5_exists_parallel(G0, p) else 6
        gam = np.log(n / k) / np.log(d)
        print(f"k={k} s={s}: n={n} dimG0={G0.shape[0]} d={'>=6' if d == 6 else d} gamma{'<=' if d == 6 else '='}{gam:.3f} [{time.time()-t0:.0f}s]", flush=True)
        out.append({"k": k, "sample": s, "n": n, "dimG0": int(G0.shape[0]), "d": d, "d_is_lower_bound": d == 6, "cols_in_X": cols})
        json.dump({"charge": "x1(x1-1)(x1-2)(x1-3)(x1-4)", "p": p, "m": m, "r": r, "seed": seed, "draws": out}, open(out_path, "w"), indent=1)
