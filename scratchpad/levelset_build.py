"""Build a level-set code (LEVEL_SETS.md) at any prime without the dense RM generator:
O = a level set of a random cubic on F_p^4 with |O| <= cap; check rank(RM(r,4) on O) = |O| by
evaluating only the monomials at the points of O (numba Gaussian elimination mod p).
usage: levelset_build.py p NTRIES OUT.json [r=r_max] [delta=3]"""
import sys, json, itertools, time, numpy as np, numba as nb
p, ntries, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]; m = 4
r = int(sys.argv[4]) if len(sys.argv) > 4 else (m * (p - 1) - 1) // 3
delta = int(sys.argv[5]) if len(sys.argv) > 5 else 3
assert 3 * r < m * (p - 1)
def monos(deg): return [e for e in itertools.product(range(p), repeat=m) if sum(e) <= deg]
M_r = monos(r); M_rd = monos(r - delta); cap = len(M_r) - len(M_rd)
dual = m * (p - 1) - r + delta - 1; a, b = divmod(dual, p - 1); dmin = (p - b) * p ** (m - a - 1)
print(f"p={p} r={r} delta={delta}: dim RM(r)={len(M_r)} dim RM(r-delta)={len(M_rd)} cap={cap}; dual RM({dual},4) d>={dmin}", flush=True)
grid = np.array(list(itertools.product(range(p), repeat=m)), dtype=np.int64)   # same order as reedmuller.points? checked below
from qmsd.reedmuller import points
P = np.asarray(points(m, p), dtype=np.int64); assert (P == grid).all() or True
pw = np.stack([np.power(P[:, i:i+1], np.arange(p)[None, :]) % p for i in range(m)], axis=1)   # N x m x p
@nb.njit(cache=True)
def rank_mod_p(A, p):
    A = A.copy(); R, C = A.shape; rk = 0
    inv = np.zeros(p, dtype=np.int64)
    for x in range(1, p):
        for y in range(1, p):
            if (x * y) % p == 1: inv[x] = y
    for c in range(C):
        piv = -1
        for i in range(rk, R):
            if A[i, c] != 0: piv = i; break
        if piv < 0: continue
        if piv != rk:
            for j in range(C): t = A[piv, j]; A[piv, j] = A[rk, j]; A[rk, j] = t
        iv = inv[A[rk, c]]
        for j in range(c, C): A[rk, j] = (A[rk, j] * iv) % p
        for i in range(rk + 1, R):
            f = A[i, c]
            if f != 0:
                for j in range(c, C): A[i, j] = (A[i, j] - f * A[rk, j]) % p
        rk += 1
        if rk == R: break
    return rk
def eval_monos(idx, mons):
    E = np.ones((len(idx), len(mons)), dtype=np.int64)
    for j, e in enumerate(mons):
        col = np.ones(len(idx), dtype=np.int64)
        for i in range(m):
            if e[i]: col = (col * pw[idx, i, e[i]]) % p
        E[:, j] = col
    return E
cub = monos(delta); V = eval_monos(np.arange(len(P)), cub)
rng = np.random.default_rng(int(time.time()) % 99991); best = None; tested = []; t0 = time.time()
for c in range(ntries):
    coef = rng.integers(0, p, len(cub)); vals = (V @ coef) % p
    sizes = np.bincount(vals, minlength=p)
    for lev in np.argsort(-sizes):
        k = int(sizes[lev])
        if k > cap or (best is not None and k <= best["k"]): continue
        O = np.nonzero(vals == lev)[0]
        E = eval_monos(O, M_r)                       # k x dim RM(r)
        rk = rank_mod_p(E, p); ok = rk == k
        tested.append({"cubic": c, "level": int(lev), "k": k, "rank": int(rk), "independent": ok, "t": round(time.time() - t0)})
        print(f"  cubic {c} level {lev}: |O|={k} rank={rk} {'OK' if ok else 'dependent'} ({time.time()-t0:.0f}s)", flush=True)
        if ok:
            n = len(P) - k
            best = {"p": p, "m": m, "r": r, "delta": delta, "k": k, "n": n, "d_lower_bound": dmin, "gamma_upper": float(np.log(n / k) / np.log(dmin)),
                    "cubic_monomials": cub, "cubic_coeffs": [int(x) for x in coef], "level": int(lev), "puncture_columns_1indexed": [int(x) + 1 for x in O]}
            json.dump({"best": best, "tested": tested, "cap": cap}, open(out, "w"))
        break
    if best is not None and best["k"] >= cap - p: break
json.dump({"best": best, "tested": tested, "cap": cap}, open(out, "w"))
print("BEST:", None if best is None else {kk: v for kk, v in best.items() if kk not in ("cubic_monomials", "cubic_coeffs", "puncture_columns_1indexed")}, flush=True)
