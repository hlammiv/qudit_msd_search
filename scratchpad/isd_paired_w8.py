"""Stern ISD for the weight-8 d''-words of a paired-cap code: weight-8 dependencies of H6 = RM(6,3) shortened at Y
(39 x 298 at |Y| = 45; no dependency of weight <= 7, cert_d2_paired.py).  Information set split 2|2, window of ELL zero
rows, 4 redundancy points.  Saturation evidence for part (ii) of census_paired_w8.py (not a proof).
Writes the distinct projective words (columns = points of F_7^3, coefficients first = 1) and the discovery curve.
usage: PYTHONPATH=. python scratchpad/isd_paired_w8.py SWEEP.json T SAMPLE ITERS SEED OUT.json [ELL=10]"""
import sys, json, time, itertools, numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.triorthogonal import shorten_matrix
from qmsd.structured_ad import _rref

src, t, s, iters, seed, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
ELL = int(sys.argv[7]) if len(sys.argv) > 7 else 10; p = 7; WT = 8
dr = [d for d in json.load(open(src))["draws"] if d["t"] == t and d["sample"] == s][0]; Y = dr["Y"]; Yset = set(Y)
offY = np.array([x for x in range(p ** 3) if x not in Yset])
R6 = np.asarray(rm_generator(6, 3, p), dtype=np.int64) % p
R, piv = _rref(np.asarray(shorten_matrix(R6, Y, p), dtype=np.int64) % p, p); H = R[:len(piv)]
r, N = H.shape; rng = np.random.default_rng(seed); found = {}; curve = []; T0 = time.time()

def pair_list(cols, Hs, first_one):
    pr = np.array(list(itertools.combinations(cols, 2)), dtype=np.int64)
    cs = [(1, b) for b in range(1, p)] if first_one else [(a, b) for a in range(1, p) for b in range(1, p)]
    V = np.concatenate([(a * Hs[:, pr[:, 0]] + b * Hs[:, pr[:, 1]]).T % p for a, b in cs])
    C = np.concatenate([np.tile(np.array([a, b]), (len(pr), 1)) for a, b in cs])
    return np.tile(pr, (len(cs), 1)), C, V

for it in range(iters):
    perm = rng.permutation(N); Rr, pv = _rref(H[:, perm], p); Hs = Rr[:r]
    if len(pv) < r: continue
    info = np.array([c for c in range(N) if c not in set(pv)]); rng.shuffle(info)
    I1, I2 = info[: len(info) // 2], info[len(info) // 2:]
    rowof = {c: i for i, c in enumerate(pv)}                       # pivot column -> its row (identity part)
    win = rng.choice(r, ELL, replace=False); rest = np.setdiff1d(np.arange(r), win)
    P1, C1, V1 = pair_list(I1, Hs, True); P2, C2, V2 = pair_list(I2, Hs, False)
    W = p ** np.arange(ELL, dtype=np.int64)
    k1 = (V1[:, win] @ W); k2 = (((-V2[:, win]) % p) @ W)
    o1 = np.argsort(k1, kind="stable"); k1s = k1[o1]
    lo = np.searchsorted(k1s, k2, "left"); hi = np.searchsorted(k1s, k2, "right")
    for j in np.nonzero(hi > lo)[0]:
        for i in o1[lo[j]:hi[j]]:
            tot = (V1[i] + V2[j]) % p                                 # = -e_R (window rows are zero)
            nzr = np.nonzero(tot[rest])[0]
            if len(nzr) != WT - 4: continue
            sup = {int(perm[P1[i, 0]]): int(C1[i, 0]), int(perm[P1[i, 1]]): int(C1[i, 1]),
                   int(perm[P2[j, 0]]): int(C2[j, 0]), int(perm[P2[j, 1]]): int(C2[j, 1])}
            for rr in rest[nzr]: sup[int(perm[pv[rr]])] = int((-tot[rr]) % p)
            cols = sorted(sup); e = np.zeros(N, dtype=np.int64); e[cols] = [sup[c] for c in cols]
            assert not ((H @ e) % p).any() and len(cols) == WT
            inv = pow(sup[cols[0]], p - 2, p)
            found[tuple((int(offY[c]), sup[c] * inv % p) for c in cols)] = it
    if it % 200 == 0 or it == iters - 1:
        curve.append((it + 1, len(found)))
        print(f"[{time.time()-T0:6.0f}s] iter {it+1}/{iters}: distinct weight-8 d''-words {len(found)}", flush=True)
        json.dump({"t": t, "sample": s, "seed": seed, "ell": ELL, "iters_done": it + 1, "curve": curve,
                   "words": [[list(c) for c in w] for w in found]}, open(out, "w"))
