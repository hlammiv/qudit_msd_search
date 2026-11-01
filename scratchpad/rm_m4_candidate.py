"""Ordinary punctured RM(r, 4) candidates at the accidental-word window edge predicted by the law
(scratchpad/window_law.py): sample k outputs with bounded plane occupancy, check full rank, record
the max outputs per plane / line / 3-flat, and save G0 columns for the MITM stages on lenore.
usage: rm_m4_candidate.py p r k seed OUT.json [max_per_plane]"""
import sys, json, random, itertools, time, numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.structured_ad import _rank
p, r, k, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); out = sys.argv[5]
cap = int(sys.argv[6]) if len(sys.argv) > 6 else None
m = 4; t0 = time.time()
P = np.asarray(points(m, p), dtype=np.int64); N = len(P)
G = np.asarray(rm_generator(r, m, p), dtype=np.int64) % p
rk = _rank(G.copy(), p); print(f"p={p} m=4 r={r}: RM rows {G.shape[0]} rank {rk}; N={N}; k={k} -> checks {rk-k}", flush=True)
# flats: all 2-dim quotient maps (planes = fibres of a rank-2 linear map), lines = fibres of rank-3 maps,
# 3-flats = fibres of rank-1 maps.  Occupancy = outputs per fibre.
def rref_rows(rows):
    M = np.array(rows) % p; piv = []
    rr = 0
    for c in range(m):
        pr = [i for i in range(rr, M.shape[0]) if M[i, c]]
        if not pr: continue
        i = pr[0]; M[[rr, i]] = M[[i, rr]]; M[rr] = (M[rr] * pow(int(M[rr, c]), -1, p)) % p
        for j in range(M.shape[0]):
            if j != rr and M[j, c]: M[j] = (M[j] - M[j, c] * M[rr]) % p
        rr += 1
        if rr == M.shape[0]: break
    return tuple(M.flatten())
normals = [v for v in itertools.product(range(p), repeat=m) if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1]
hyper = [np.array(v) for v in normals]
planes = {rref_rows([normals[i], normals[j]]) for i in range(len(normals)) for j in range(i + 1, len(normals))}
planes = [np.array(s).reshape(2, m) for s in planes]
print(f"{len(hyper)} hyperplane directions, {len(planes)} plane-quotient maps ({time.time()-t0:.0f}s)", flush=True)
def max_occupancy(mask, maps):
    best = 0
    for S in maps:
        key = (P @ S.T) % p
        key = key[:, 0] if key.ndim == 1 or key.shape[1] == 1 else sum(key[:, i] * p ** i for i in range(key.shape[1]))
        cnt = np.bincount(key[mask], minlength=p ** S.shape[0] if S.ndim == 2 else p)
        best = max(best, int(cnt.max()))
    return best
rng = random.Random(seed)
for attempt in range(50):
    cols = sorted(rng.sample(range(N), k)); mask = np.zeros(N, dtype=bool); mask[cols] = True
    occ_plane = max_occupancy(mask, planes)
    if cap is not None and occ_plane > cap: continue
    if _rank(G[:, cols].T.copy(), p) < k: continue
    break
else:
    raise SystemExit("no sample met the constraints")
occ_hyper = max_occupancy(mask, [h.reshape(1, m) for h in hyper])
dimA = rk - k
json.dump({"p": p, "m": m, "r": r, "k": k, "n": N - k, "seed": seed, "attempts": attempt + 1, "rank": rk, "dim_A": dimA,
           "max_outputs_per_plane": occ_plane, "mean_outputs_per_plane": k * 1.0 / p ** 2, "max_outputs_per_hyperplane": occ_hyper,
           "gamma_if_d": {d: float(np.log((N - k) / k) / np.log(d)) for d in (5, 6, 7, 8)}, "puncture_columns_1indexed": [c + 1 for c in cols],
           "seconds": round(time.time() - t0)}, open(out, "w"), indent=1)
print(f"saved {out}: n={N-k} dim A={dimA} max/plane {occ_plane} (mean {k/p**2:.1f}) max/hyperplane {occ_hyper} gamma(d=6)={np.log((N-k)/k)/np.log(6):.3f} ({time.time()-t0:.0f}s)", flush=True)
