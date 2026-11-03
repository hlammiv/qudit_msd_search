"""Level-set codes at large p (29, 31): build the largest cubic level set under the cap and certify that its
RM(r,4) columns are independent with a parallel, memory-lean elimination (int16 storage, exact modular powers).
Two independent routines: --method row (row echelon, columns left to right) and --method col (transposed,
columns right to left).  usage: levelset_big.py p delta OUT.json --method row|col [--from CODE.json]"""
import sys, json, itertools, time, numpy as np, numba as nb
from qmsd.reedmuller import points
p, delta, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]; m = 4
method = sys.argv[sys.argv.index("--method") + 1]; src = sys.argv[sys.argv.index("--from") + 1] if "--from" in sys.argv else None
r = (m * (p - 1) - 1) // 3; t0 = time.time()
POW = np.array([[pow(x, k, p) for k in range(p)] for x in range(p)], dtype=np.int64)
P = np.asarray(points(m, p), dtype=np.int64)
def monos(deg): return [e for e in itertools.product(range(p), repeat=m) if sum(e) <= deg]
M_r = monos(r); cap = len(M_r) - len(monos(r - delta))
dual = m * (p - 1) - r + delta - 1; a, b = divmod(dual, p - 1); dmin = (p - b) * p ** (m - a - 1)
if src:
    bb = json.load(open(src))["best"]; O = np.array(sorted(c - 1 for c in bb["puncture_columns_1indexed"])); coef = bb["cubic_coeffs"]; lev = bb["level"]; cub = [tuple(e) for e in bb["cubic_monomials"]]
else:
    cub = monos(delta); rng = np.random.default_rng(int(time.time()) % 99991)
    V = np.ones((len(P), len(cub)), dtype=np.int64)
    for j, e in enumerate(cub):
        for i in range(m):
            if e[i]: V[:, j] = (V[:, j] * POW[P[:, i], e[i]]) % p
    while True:
        coef = rng.integers(0, p, len(cub)); vals = (V @ coef) % p; sizes = np.bincount(vals, minlength=p)
        ok = [l for l in range(p) if sizes[l] <= cap]
        if ok: lev = max(ok, key=lambda l: sizes[l]); O = np.nonzero(vals == lev)[0]; break
    coef = [int(x) for x in coef]
k = len(O); n = len(P) - k
print(f"p={p} r={r} delta={delta}: cap {cap}, |O| = {k}, dim RM(r) = {len(M_r)}; dual RM({dual},4) d >= {dmin}; gamma <= {np.log(n/k)/np.log(dmin):.4f} ({time.time()-t0:.0f}s)", flush=True)
# level-set check (recompute from the stored polynomial)
vals = np.zeros(len(P), dtype=np.int64)
for e, c in zip(cub, coef):
    if c:
        t = np.full(len(P), c, dtype=np.int64)
        for i in range(m): t = (t * POW[P[:, i], e[i]]) % p
        vals = (vals + t) % p
if not src:
    json.dump({"best": {"p": p, "m": m, "r": r, "delta": delta, "k": int(k), "n": int(n), "d_lower_bound": int(dmin), "gamma_upper": float(np.log(n / k) / np.log(dmin)),
                        "cubic_monomials": [list(e) for e in cub], "cubic_coeffs": [int(x) for x in coef], "level": int(lev), "puncture_columns_1indexed": [int(x) + 1 for x in O]}},
              open(out.replace(".json", "_set.json"), "w"))
print(f"level set recomputed equals O: {(np.nonzero(vals == lev)[0] == O).all() and len(np.nonzero(vals == lev)[0]) == k}", flush=True)
E = np.ones((k, len(M_r)), dtype=np.int16)
for j, e in enumerate(M_r):
    col = np.ones(k, dtype=np.int64)
    for i in range(m):
        if e[i]: col = (col * POW[P[O, i], e[i]]) % p
    E[:, j] = col.astype(np.int16)
print(f"E built {E.shape} ({E.nbytes/1e9:.1f} GB, {time.time()-t0:.0f}s)", flush=True)

@nb.njit(parallel=True, cache=True)
def rank_rows(A, p, every):
    R, C = A.shape; inv = np.zeros(p, dtype=np.int64)
    for x in range(1, p):
        for y in range(1, p):
            if (x * y) % p == 1: inv[x] = y
    rk = 0
    for c in range(C):
        piv = -1
        for i in range(rk, R):
            if A[i, c] != 0: piv = i; break
        if piv < 0: continue
        if piv != rk:
            for j in range(C): t = A[piv, j]; A[piv, j] = A[rk, j]; A[rk, j] = t
        iv = inv[A[rk, c]]
        for j in range(c, C): A[rk, j] = (A[rk, j] * iv) % p
        for i in nb.prange(rk + 1, R):
            f = np.int64(A[i, c])
            if f != 0:
                for j in range(c, C): A[i, j] = (np.int64(A[i, j]) - f * np.int64(A[rk, j])) % p
        rk += 1
        if every > 0 and rk % every == 0: print('  pivots', rk, 'of', R)
        if rk == R: break
    return rk

@nb.njit(parallel=True, cache=True)
def rank_cols_rev(B, p, every):
    # B = E^T (C x R): eliminate on the transpose, pivot columns right to left
    R, C = B.shape; inv = np.zeros(p, dtype=np.int64)
    for x in range(1, p):
        for y in range(1, p):
            if (x * y) % p == 1: inv[x] = y
    rk = 0
    for c in range(C - 1, -1, -1):
        piv = -1
        for i in range(rk, R):
            if B[i, c] != 0: piv = i; break
        if piv < 0: continue
        if piv != rk:
            for j in range(C): t = B[piv, j]; B[piv, j] = B[rk, j]; B[rk, j] = t
        iv = inv[B[rk, c]]
        for i in nb.prange(rk + 1, R):
            f = (np.int64(B[i, c]) * iv) % p
            if f != 0:
                for j in range(C): B[i, j] = (np.int64(B[i, j]) - f * np.int64(B[rk, j])) % p
        rk += 1
        if every > 0 and rk % every == 0: print('  pivots', rk, 'of', R)
        if rk == R: break
    return rk

EXTRA = 40
seed = 1 if method == "row" else 2                     # different column subsets for the two routines
cols = np.sort(np.random.default_rng(seed).choice(E.shape[1], size=min(E.shape[1], int(1.25 * k) + EXTRA), replace=False))
Es = np.ascontiguousarray(E[:, cols]); print(f"column-subset certificate: {Es.shape[1]} of {E.shape[1]} columns (seed {seed}); full rank on the subset implies full rank of E", flush=True)
import threading
t_start = time.time()
if method == "row": rk = rank_rows(Es, p, 1000)
else: rk = rank_cols_rev(np.ascontiguousarray(Es.T), p, 1000)
if rk < k:
    print(f"subset rank {rk} < {k}: falling back to all columns", flush=True)
    if method == "row": rk = rank_rows(E, p, 1000)
    else: rk = rank_cols_rev(np.ascontiguousarray(E.T), p, 1000)
print(f"rank ({method}): {rk} of {k} -> {'INDEPENDENT' if rk == k else 'DEPENDENT'} ({time.time()-t0:.0f}s)", flush=True)
json.dump({"best": {"p": p, "m": m, "r": r, "delta": delta, "k": int(k), "n": int(n), "d_lower_bound": int(dmin), "gamma_upper": float(np.log(n / k) / np.log(dmin)),
                    "cubic_monomials": [list(e) for e in cub], "cubic_coeffs": [int(x) for x in coef], "level": int(lev), "puncture_columns_1indexed": [int(x) + 1 for x in O],
                    "rank": int(rk), "method": method}}, open(out, "w"))
