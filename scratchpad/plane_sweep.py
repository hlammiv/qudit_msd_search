"""Multi-seed sweep, cheap stages only (fits a 2 GB box): for each random cube-class output set
of size k on the p=7 m=4 r=7 sphere charge, compute the plane-restricted distance d_plane
(rigorous upper bound on d; exact on 8 of 9 known draws) and, if d_plane >= 7, the class-block
weight <= 5 check. Survivors (d_plane >= 7 and no block word of weight <= 5) are candidates for
the full certificate on lenore. Appends JSON lines; resumable.
usage: plane_sweep.py K SEED0 NSEEDS OUT.jsonl"""
import sys, json, random, itertools, time, os, numpy as np
from qmsd.reedmuller import points
from qmsd.weighted import weighted_support, quadric, class_points
from qmsd.structured_ad import _rank, _null_space
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import weight_d_exists_fast
k, seed0, nseeds, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
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
COEFS = {nu: np.array(list(itertools.product(range(p), repeat=nu)))[1:] for nu in range(1, 7)}
def min_null_weight(M, cap):
    ns = np.asarray(_null_space(M, p)) % p; nu = ns.shape[0]
    if nu == 0: return None
    if nu <= 6:
        return int(((COEFS[nu] @ ns) % p != 0).sum(axis=1).min())
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
    best = 99
    for Sb in subs:
        v = (PX @ Sb.T) % p; key = v[:, 0] * p + v[:, 1]
        order = np.argsort(key, kind="stable"); ks = key[order]
        starts = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]]); ends = np.r_[starts[1:], len(ks)]
        for a, b in zip(starts, ends):
            idx = order[a:b]; cc = local[idx]; cc = cc[cc >= 0]
            if len(cc) < 2: continue
            M = Gc[:, cc]
            if _rank(M.T.copy(), p) == len(cc): continue
            w = min_null_weight(M, best)
            if w is not None and w < best:
                best = w
                if best < 7: return best, G0, Gc
    return best, G0, Gc
done = set()
if os.path.exists(out):
    for line in open(out):
        if line.strip(): rec = json.loads(line); done.add((rec["k"], rec["seed"]))
for seed in range(seed0, seed0 + nseeds):
    if (k, seed) in done: continue
    rng = random.Random(1000 * k + seed)
    cols = sorted(rng.sample(cand, k))
    while _rank(GX[:, cols].T.copy(), p) < k: cols = sorted(rng.sample(cand, k))
    t0 = time.time(); dpl, G0, Gc = plane_distance(cols)
    rec = {"k": k, "seed": seed, "n": G0.shape[1], "d_plane": dpl, "t_plane": round(time.time() - t0)}
    if dpl >= 7:
        t1 = time.time(); db = None
        for dd in (2, 3, 4, 5):
            if weight_d_exists_fast(Gc, p, dd): db = dd; break
        rec["block_d_le5"] = db; rec["t_block"] = round(time.time() - t1); rec["cols_in_X"] = cols
    with open(out, "a") as f: f.write(json.dumps(rec) + "\n")
    print(json.dumps({kk: vv for kk, vv in rec.items() if kk != "cols_in_X"}), flush=True)
