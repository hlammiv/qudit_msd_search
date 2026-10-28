"""Exact flat cap for a weighted support: draw random class output sets of size k and report the
minimum weight of a code word supported on an affine codim-1 flat (m = 2 lines, m = 3 planes) or
codim-2 flat (m = 4 planes), using MITM to weight 6 and the 3+4 witness for weight 7 (small flats
only). Words on a flat are code words, so this is a rigorous upper bound on d.
usage: p m r CHARGE UNIT K NSEEDS   CHARGE = q:1,0:+x | q:1,1:+0 | cubic:+0 | lin:0,0,1:+0"""
import sys, random, itertools, time, numpy as np
from qmsd.reedmuller import points
from qmsd.weighted import weighted_support, quadric, class_points
from qmsd.structured_ad import _rank
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import weight_d_exists_fast, weight7_witness_3x4
p, m, r = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); spec, unit, k, nseeds = sys.argv[4], int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
P = np.asarray(points(m, p), dtype=np.int64)
kind, *rest = spec.split(":")
if kind == "q":
    diag = tuple(int(x) for x in rest[0].split(",")); tail = rest[1]
    g = quadric(p, m, diag, lin=tuple(1 if d == 0 else 0 for d in diag)) if tail == "+x" else quadric(p, m, diag, const=int(tail))
elif kind == "cubic":
    g = (np.sum(P ** 3, axis=1) + int(rest[0])) % p
elif kind == "lin":
    lin = np.array([int(x) for x in rest[0].split(",")]); g = (P @ lin + int(rest[1])) % p
GX, X, mu = weighted_support(p, m, g, r); cand = class_points(mu, p, unit)
print(f"p={p} m={m} r={r} charge {spec} unit {unit}: |X|={len(X)} rank={_rank(GX, p)} N_c={len(cand)} rho={_rank(GX[:, cand], p)}", flush=True)
PX = P[X]
normals = [v for v in itertools.product(range(p), repeat=m) if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1]
codim = int(sys.argv[8]) if len(sys.argv) > 8 else (1 if m <= 3 else 2)
def flats():
    if codim == 1:
        for a in normals:
            key = (PX @ np.array(a)) % p
            for c in range(p): yield np.flatnonzero(key == c)
    elif codim == 2:
        for a, b in itertools.combinations(normals, 2):
            M = np.array([a, b]) % p
            if _rank(M.copy(), p) < 2: continue
            v = (PX @ M.T) % p; key = v[:, 0] * p + v[:, 1]
            order = np.argsort(key, kind="stable"); ks = key[order]
            starts = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]]); ends = np.r_[starts[1:], len(ks)]
            for s, e in zip(starts, ends): yield order[s:e]
    else:  # lines (codim m-1): direction v, group by the projection orthogonal to v
        for v in normals:
            Mv = np.array([w for w in normals if _rank(np.array([v, w]) % p, p) == 2][:m - 1]) % p
            # complete to a basis of normals orthogonal... simpler: group by all m-1 independent normals vanishing on v
            comp = [w for w in itertools.product(range(p), repeat=m) if any(w) and (np.array(w) @ np.array(v)) % p == 0]
            B = []
            for w in comp:
                if _rank(np.array(B + [w]) % p, p) == len(B) + 1: B.append(w)
                if len(B) == m - 1: break
            key = np.zeros(len(PX), dtype=np.int64)
            for w in B: key = key * p + (PX @ np.array(w)) % p
            order = np.argsort(key, kind="stable"); ks = key[order]
            starts = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]]); ends = np.r_[starts[1:], len(ks)]
            for s, e in zip(starts, ends): yield order[s:e]
rng = random.Random(2026)
for s in range(nseeds):
    cols = sorted(rng.sample(cand, k))
    while _rank(GX[:, cols].T.copy(), p) < k: cols = sorted(rng.sample(cand, k))
    G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p
    keep = [l for l in range(GX.shape[1]) if l not in set(cols)]; pos = {l: j for j, l in enumerate(keep)}
    t0 = time.time(); best = 99; nfl = 0
    for idx in flats():
        cc = np.array([pos[l] for l in idx if l in pos], dtype=np.int64)
        if len(cc) < 2: continue
        M = np.ascontiguousarray(G0[:, cc]); M = M[[t for t in range(M.shape[0]) if M[t].any()]]
        if M.shape[0] == 0: best = 1; break
        if _rank(M.T.copy(), p) == len(cc): continue
        nfl += 1
        for d in range(2, min(best, 7)):
            if weight_d_exists_fast(M, p, d): best = d; break
        else:
            if best > 7 and len(cc) <= 60 and weight7_witness_3x4(M, p) is not None: best = 7
        if best <= 5: break
    print(f"seed {s}: k={k} n={G0.shape[1]} checks={G0.shape[0]} d_flat={'<=' if best == 99 else ''}{min(best, 99) if best < 99 else '>=8'} dependent flats {nfl} [{time.time()-t0:.0f}s]", flush=True)
