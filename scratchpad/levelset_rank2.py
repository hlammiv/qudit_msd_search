"""Second rank routine for level-set codes (corrected exact powers): column-major elimination written
independently of fast_matrices._rref and levelset_build.rank_mod_p.  usage: levelset_rank2.py CODE.json"""
import sys, json, itertools, time, numpy as np, numba as nb
from qmsd.reedmuller import points
@nb.njit(cache=True)
def rank_colmajor(A, p):
    # eliminate on the transpose, pivoting from the LAST column: different order and data layout
    M = A.T.copy(); R, C = M.shape; inv = np.zeros(p, dtype=np.int64)
    for x in range(1, p):
        for y in range(1, p):
            if (x * y) % p == 1: inv[x] = y
    rk = 0
    for c in range(C - 1, -1, -1):
        piv = -1
        for i in range(rk, R):
            if M[i, c] % p != 0: piv = i; break
        if piv < 0: continue
        if piv != rk:
            for j in range(C): t = M[piv, j]; M[piv, j] = M[rk, j]; M[rk, j] = t
        iv = inv[M[rk, c] % p]
        for i in range(rk + 1, R):
            f = (M[i, c] * iv) % p
            if f != 0:
                for j in range(C): M[i, j] = (M[i, j] - f * M[rk, j]) % p
        rk += 1
        if rk == R: break
    return rk
b = json.load(open(sys.argv[1]))["best"]; p, m, r = b["p"], b["m"], b["r"]; t0 = time.time()
POW = np.array([[pow(x, k, p) for k in range(p)] for x in range(p)], dtype=np.int64)
P = np.asarray(points(m, p), dtype=np.int64); O = np.array(sorted(c - 1 for c in b["puncture_columns_1indexed"]))
mons = [e for e in itertools.product(range(p), repeat=m) if sum(e) <= r]
E = np.ones((len(O), len(mons)), dtype=np.int64)
for j, e in enumerate(mons):
    for i in range(m):
        if e[i]: E[:, j] = (E[:, j] * POW[P[O, i], e[i]]) % p
rk = rank_colmajor(np.ascontiguousarray(E), p)
print(f"rank (column-major, exact powers): {rk} (|O| = {len(O)}) -> {'INDEPENDENT' if rk == len(O) else 'DEPENDENT'} ({time.time()-t0:.0f}s)", flush=True)
