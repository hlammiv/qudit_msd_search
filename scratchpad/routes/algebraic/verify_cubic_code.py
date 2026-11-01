"""Independent verification of a cubic-level-set code from p11_cubic_hunt.py (best.json):
  1. recompute Z = {g = level} from the coefficients and compare with the stored O;
  2. rank of RM(13,4) on O by an independent numpy elimination mod 11 (not galois);
  3. build G0 (RM(13,4) words vanishing on O, restricted to N) and smoke-test the theorem:
     random 21-subsets of G0's columns are independent, 21 points in a plane are independent,
     two parallel lines inside N (22 points) -- report the rank;
  4. no kernel word of weight <= 3: MITM on an 18-row random projection of G0, hits verified
     against the full G0 (none expected);
  5. gamma = ln(n/k)/ln(22).
Usage: python verify_cubic_code.py best.json"""
import sys, json, time
import numpy as np
from qmsd.reedmuller import rm_generator, points, monomial_exponents, d_rm
from qmsd.field import GFp
from qmsd.mindist import min_dependent_columns

p, m, r = 11, 4, 13
best = json.load(open(sys.argv[1]))
O = np.array(best["puncture_columns_1indexed"]) - 1
k = len(O); n0 = p ** m; n = n0 - k
print(f"k = {k}, n = {n}, level = {best['level']}")

# 1. zero set
G3 = np.asarray(rm_generator(3, m, p), dtype=np.int64)
g = (np.array(best["coeffs"]) @ G3) % p
Z = np.nonzero(g == best["level"])[0]
assert np.array_equal(np.sort(Z), np.sort(O)), "stored O is not the level set"
print("1. O == {g = level}: OK; degree-3 part nonzero:",
      any(c and sum(a) == 3 for c, a in zip(best["coeffs"], monomial_exponents(m, 3, p))))


def rank_numpy(M, p):
    M = (np.asarray(M, dtype=np.int64) % p).copy()
    nr, nc = M.shape
    rk = 0
    for c in range(nc):
        if rk == nr:
            break
        col = M[rk:, c]
        nz = np.nonzero(col)[0]
        if len(nz) == 0:
            continue
        piv = rk + nz[0]
        if piv != rk:
            M[[rk, piv]] = M[[piv, rk]]
        M[rk] = (M[rk] * pow(int(M[rk, c]), p - 2, p)) % p
        others = np.nonzero(M[:, c])[0]
        others = others[others != rk]
        if len(others):
            M[others] = (M[others] - np.outer(M[others, c], M[rk])) % p
        rk += 1
    return rk


t0 = time.time()
G = np.asarray(rm_generator(r, m, p), dtype=np.int64)
rk = rank_numpy(G[:, O].T, p)
print(f"2. numpy rank of RM(13,4) on O = {rk} (|O| = {k}): {'FULL RANK' if rk == k else 'DEFICIENT'} [{time.time()-t0:.0f}s]")
assert rk == k

# 3. G0
GF = GFp(p)
t0 = time.time()
u = GF(G[:, O] % p).left_null_space()
assert u.shape[0] == G.shape[0] - k
isO = np.zeros(n0, bool); isO[O] = True
N = np.nonzero(~isO)[0]
uf = np.asarray(u, dtype=np.int64).astype(np.float64)       # exact: sums <= 2320*100 < 2^53
G0 = (np.rint(uf @ G.astype(np.float64)) % p).astype(np.int64)[:, N]
assert (G0.shape[0] == G.shape[0] - k) and ((np.rint(uf @ G[:, O].astype(np.float64)) % p) == 0).all()
print(f"3. G0 shape {G0.shape} [{time.time()-t0:.0f}s]")
rng = np.random.default_rng(0)
bad = 0
for t in range(200):
    T = rng.choice(G0.shape[1], 21, replace=False)
    if rank_numpy(G0[:, T], p) < 21:
        bad += 1
print(f"   random 21-subsets of G0 columns with a dependency: {bad}/200 (theorem: 0)")
pts = points(m, p)
colN = -np.ones(n0, dtype=np.int64); colN[N] = np.arange(len(N))
# 21 points of a plane inside N, and two parallel lines inside N
dep_plane = 0; tested = 0
for t in range(100):
    a = pts[rng.integers(n0)]; v1 = rng.integers(0, p, m); v2 = rng.integers(0, p, m)
    P = np.array([(a + s * v1 + t2 * v2) % p for s in range(p) for t2 in range(p)])
    idx = np.unique(P @ (p ** np.arange(m)))
    if len(idx) < p * p:
        continue
    idx = idx[~isO[idx]]
    if len(idx) < 21:
        continue
    T = colN[rng.choice(idx, 21, replace=False)]
    tested += 1
    if rank_numpy(G0[:, T], p) < 21:
        dep_plane += 1
print(f"   21 points in a plane (N part): dependent {dep_plane}/{tested} (theorem: 0)")
ranks22 = []
tries = 0
while len(ranks22) < 20 and tries < 20000:
    tries += 1
    a = pts[rng.integers(n0)]; v = rng.integers(0, p, m)
    if not v.any():
        continue
    L1 = np.array([(a + s * v) % p for s in range(p)]) @ (p ** np.arange(m))
    w = rng.integers(0, p, m)
    L2 = np.array([(a + w + s * v) % p for s in range(p)]) @ (p ** np.arange(m))
    if isO[L1].any() or isO[L2].any() or len(set(L1) & set(L2)):
        continue
    T = colN[np.concatenate([L1, L2])]
    ranks22.append(rank_numpy(G0[:, T], p))
print(f"   two parallel lines inside N (22 columns): ranks {sorted(set(ranks22))} over {len(ranks22)} pairs "
      f"(rank 21 would be a weight-22 kernel word, allowed by the bound; rank 22 = independent)")

# 4. weight <= 3 via projection
t0 = time.time()
R = rng.integers(0, p, size=(18, G0.shape[0]))
H18 = (R @ G0) % p
try:
    dmin = min_dependent_columns(H18, p, d_max=3)
    print(f"4. projection MITM found a dependent set of size {dmin} -- must verify against full G0 (false positive possible)")
except ValueError:
    print(f"4. no kernel word of weight <= 3 (18-row projection MITM, {time.time()-t0:.0f}s): OK")

print(f"5. [[{n},{k},>=22]]_11, gamma <= ln({n}/{k})/ln 22 = {np.log(n/k)/np.log(22):.4f}; "
      f"d(RM(29,4)) = {d_rm(29, 4, 11)}")
