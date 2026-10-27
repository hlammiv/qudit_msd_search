"""Plane-restricted distance: for every affine plane, the surviving class columns of G0 restricted
to it; if they are dependent, the min weight of the null space (nullity is tiny: <= 21 columns
vs 21-31 checks).  Words supported inside a plane are words of the code, so this is an upper
bound on d; validated against the certified distances.  3-flats spot-checked with numba MITM."""
import json, random, itertools, time, sys, numpy as np
from qmsd.reedmuller import points
from qmsd.weighted import weighted_support, quadric, class_points
from qmsd.structured_ad import _rank, _null_space
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import weight_d_exists_fast
p, m, r = 7, 4, 7
P = np.asarray(points(m, p), dtype=np.int64)
g = quadric(p, m, (1, 1, 1, 1)); GX, X, mu = weighted_support(p, m, g, r); cand = class_points(mu, p, 1)
normals = np.array([v for v in itertools.product(range(p), repeat=m) if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1])
def rref2(a, b):
    M = np.array([a, b]) % p
    for c in range(m):
        piv = [i for i in range(2) if M[i, c]]
        if not piv: continue
        i = piv[0]
        if i: M[[0, i]] = M[[i, 0]]
        M[0] = (M[0] * pow(int(M[0, c]), -1, p)) % p; M[1] = (M[1] - M[1, c] * M[0]) % p; break
    for c in range(m):
        if M[1, c]:
            M[1] = (M[1] * pow(int(M[1, c]), -1, p)) % p; M[0] = (M[0] - M[0, c] * M[1]) % p; break
    return tuple(M.flatten())
subs = [np.array(s).reshape(2, m) for s in {rref2(normals[i], normals[j]) for i in range(len(normals)) for j in range(i + 1, len(normals))}]
Xc = np.array([int(X[l]) for l in cand]); PX = P[Xc]
def min_null_weight(M, cap):
    """Exact min weight of the null space of M when its nullity is <= 6 (7^nu - 1 words,
    vectorised); otherwise MITM for weights 2..6 and 7 meaning '>= 7' (never happens on the
    draws with >= 16 checks)."""
    ns = np.asarray(_null_space(M, p)) % p
    nu = ns.shape[0]
    if nu == 0: return None
    if nu <= 6:
        coefs = np.array(list(itertools.product(range(p), repeat=nu)))[1:]
        W = (coefs @ ns) % p
        return int((W != 0).sum(axis=1).min())
    for dd in range(2, min(cap, 7)):
        if weight_d_exists_fast(M, p, dd): return dd
    return 7
def plane_distance(cols):
    G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p
    keep = [l for l in range(GX.shape[1]) if l not in set(cols)]; pos = {l: j for j, l in enumerate(keep)}
    sel = set(cols); surv = np.array([cand[i] not in sel for i in range(len(cand))])
    colidx = np.array([pos.get(cand[i], -1) for i in range(len(cand))])
    cls_cols = colidx[surv]; Gc = G0[:, cls_cols]; Gc = Gc[[t for t in range(Gc.shape[0]) if Gc[t].any()]]
    local = -np.ones(len(cand), dtype=int); local[np.nonzero(surv)[0]] = np.arange(surv.sum())
    best, nplanes, ndep = 99, 0, 0
    for Sb in subs:
        v = (PX @ Sb.T) % p; key = v[:, 0] * p + v[:, 1]
        order = np.argsort(key, kind="stable"); ks = key[order]
        starts = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]]); ends = np.r_[starts[1:], len(ks)]
        for a, b in zip(starts, ends):
            idx = order[a:b]; cc = local[idx]; cc = cc[cc >= 0]
            if len(cc) < 2: continue
            nplanes += 1
            M = Gc[:, cc]
            if _rank(M.T.copy(), p) == len(cc): continue
            ndep += 1
            w = min_null_weight(M, best)
            if w is not None and w < best: best = w
    return best, nplanes, ndep, G0
def sample(seed, klist):
    rng = random.Random(seed); out = []
    for k in klist:
        cols = sorted(rng.sample(cand, k))
        while _rank(GX[:, cols].T.copy(), p) < k: cols = sorted(rng.sample(cand, k))
        out.append(cols)
    return out
s2026 = sample(2026, [290, 280]); s7 = sample(7, [275, 275, 265, 265, 260, 260])
sets = [("k=290 d=2", s2026[0], 2), ("k=280 d=4", s2026[1], 4), ("k=275 s0 d_block=4", s7[0], 4), ("k=275 s1 d_block=5", s7[1], 5),
        ("k=265 s0 d_block=3", s7[2], 3), ("k=265 s1 d=6", s7[3], 6), ("k=260 s0 d_block=4", s7[4], 4), ("k=260 s1 d=7", s7[5], 7),
        ("k=270 d=6", json.load(open("runs/weighted/p7_m4_sphere_1746_270_cols.json"))["cols_in_X"], 6)]
sets = sets[-1:] + sets[7:8] + sets[5:6] + sets[:5] + sets[6:7]
res = []
for name, cols, dtrue in sets:
    t = time.time(); dpl, npl, ndep, G0 = plane_distance(cols)
    print(f"{name:22s} plane-distance {dpl:2d}  (true {dtrue})  [{npl} planes, {ndep} dependent, {time.time()-t:.0f}s]", flush=True)
    res.append({"set": name, "d_true": dtrue, "d_plane": dpl})
json.dump(res, open("runs/weighted/p7_m4_sphere_plane_distance_check.json", "w"), indent=1)
