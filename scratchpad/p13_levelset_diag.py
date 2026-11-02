"""p = 13 level-set member: search structured cubics with few zeros (diagonal a.x^3 = c, plus diagonal
with one cross term) for a level set of size <= cap = 1996 whose RM(15,4) columns are independent.
Exact powers throughout.  usage: p13_levelset_diag.py OUT.json"""
import sys, json, itertools, time, numpy as np, numba as nb
from qmsd.reedmuller import points
p, m, r, delta = 13, 4, 15, 3
POW = np.array([[pow(x, k, p) for k in range(p)] for x in range(p)], dtype=np.int64)
P = np.asarray(points(m, p), dtype=np.int64)
def monos(deg): return [e for e in itertools.product(range(p), repeat=m) if sum(e) <= deg]
M_r, M_rd = monos(r), monos(r - delta); cap = len(M_r) - len(M_rd)
@nb.njit(cache=True)
def rank_mod_p(A, p):
    A = A.copy(); R, C = A.shape; inv = np.zeros(p, dtype=np.int64)
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
        for i in range(rk + 1, R):
            f = A[i, c]
            if f != 0:
                for j in range(c, C): A[i, j] = (A[i, j] - f * A[rk, j]) % p
        rk += 1
        if rk == R: break
    return rk
def E_of(O):
    E = np.ones((len(O), len(M_r)), dtype=np.int64)
    for j, e in enumerate(M_r):
        for i in range(m):
            if e[i]: E[:, j] = (E[:, j] * POW[P[O, i], e[i]]) % p
    return E
cubes = (P ** 3) % p
t0 = time.time(); seen = set(); tested = []; best = None
# coefficient classes: cube classes mod 13 are {1}, {2}, {4} cosets of cubes -> representatives 1, 2, 4
reps = [1, 2, 4]
cands = []
for a in itertools.product(reps, repeat=4):
    vals = (cubes @ np.array(a)) % p
    for c in range(p):
        k = int((vals == c).sum())
        cands.append((k, a, c, None))
for a in itertools.product(reps, repeat=3):     # x1^3 + x2^3 + x3^3 type plus a cross term x1 x2 x4
    for cr in (1, 2):
        vals = (cubes[:, :3] @ np.array(a) + cr * P[:, 0] * P[:, 1] * P[:, 3]) % p
        for c in range(p):
            cands.append((int((vals == c).sum()), a, c, cr))
sizes = sorted({k for k, *_ in cands})
print(f"cap {cap}; {len(cands)} structured level sets; sizes {sizes[:12]} ... {sizes[-5:]}; <= cap: {sum(1 for k,*_ in cands if k <= cap)}", flush=True)
for k, a, c, cr in sorted(cands, key=lambda x: -x[0]):
    if k > cap or k in seen and len(tested) > 40: continue
    vals = ((cubes @ np.array(a)) % p) if cr is None else ((cubes[:, :3] @ np.array(a) + cr * P[:, 0] * P[:, 1] * P[:, 3]) % p)
    O = np.nonzero(vals == c)[0]
    rk = rank_mod_p(E_of(O), p); ok = rk == k; seen.add(k)
    tested.append({"k": k, "coeffs": list(a), "cross": cr, "level": c, "rank": int(rk), "independent": ok})
    print(f"  {a} cross={cr} level {c}: |O|={k} rank={rk} {'OK' if ok else 'dependent'} ({time.time()-t0:.0f}s)", flush=True)
    if ok:
        n = p ** m - k; dual = m * (p - 1) - r + delta - 1; aa, bb = divmod(dual, p - 1); d = (p - bb) * p ** (m - aa - 1)
        best = {"p": p, "m": m, "r": r, "delta": delta, "k": k, "n": n, "d_lower_bound": d, "gamma_upper": float(np.log(n / k) / np.log(d)),
                "form": "diagonal cubic" if cr is None else "diagonal + cross term x1 x2 x4", "coeffs": list(a), "cross": cr, "level": int(c),
                "puncture_columns_1indexed": [int(x) + 1 for x in O]}
        break
    if len(tested) >= 60: break
json.dump({"best": best, "tested": tested, "cap": cap}, open(sys.argv[1], "w"), indent=1)
print("BEST:", None if best is None else {kk: v for kk, v in best.items() if kk != "puncture_columns_1indexed"}, flush=True)
