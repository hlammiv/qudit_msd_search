"""Hunt for a cubic zero set in F_p^4 of size <= cap with independent RM(r,4) columns (the level-set
construction of LEVEL_SETS.md).  usage: levelset_hunt.py p r NCUBICS OUT.json [delta=3]"""
import sys, json, itertools, time, numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.structured_ad import _rank
p, r, ncub, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]; delta = int(sys.argv[5]) if len(sys.argv) > 5 else 3
m = 4; t0 = time.time()
P = np.asarray(points(m, p), dtype=np.int64); N = len(P)
G = np.asarray(rm_generator(r, m, p), dtype=np.int64) % p
Gd = np.asarray(rm_generator(r - delta, m, p), dtype=np.int64) % p
cap = G.shape[0] - Gd.shape[0]
dual = m * (p - 1) - r + delta - 1; a, b = divmod(dual, p - 1); dmin = (p - b) * p ** (m - a - 1)
print(f"p={p} r={r} delta={delta}: dim RM(r)={G.shape[0]} dim RM(r-delta)={Gd.shape[0]} cap={cap}; dual RM({dual},4) d>={dmin} ({time.time()-t0:.0f}s)", flush=True)
mons = [e for e in itertools.product(range(delta + 1), repeat=m) if sum(e) <= delta]
V = np.zeros((N, len(mons)), dtype=np.int64)
for j, e in enumerate(mons):
    col = np.ones(N, dtype=np.int64)
    for i in range(m): col = (col * (np.power(P[:, i], e[i]) % p)) % p
    V[:, j] = col
rng = np.random.default_rng(int(time.time()) % 100000); best = None; tested = []
for c in range(ncub):
    # alternate dense and sparse cubics: sparse / structured ones have more extreme zero counts
    coef = rng.integers(0, p, len(mons)); coef[0] = 0
    if c % 3:
        s = int(rng.integers(3, 9)); keep = rng.choice(len(mons) - 1, size=s, replace=False) + 1
        mask = np.zeros(len(mons), dtype=bool); mask[keep] = True; coef = np.where(mask, coef, 0)
        if not any(coef[j] and sum(mons[j]) == 3 for j in range(len(mons))): coef[len(mons) - 1] = 1 + rng.integers(0, p - 1)
    vals = (V @ coef) % p
    for lev in range(p):
        O = np.nonzero(vals == lev)[0]; k = len(O)
        if k > cap or k < cap - 3 * p: continue
        rk = _rank(G[:, O].T.copy(), p); ok = rk == k
        tested.append({"cubic": c, "level": lev, "k": k, "rank": rk, "independent": ok})
        print(f"  cubic {c} level {lev}: |O|={k} rank={rk} {'OK' if ok else 'dependent'} ({time.time()-t0:.0f}s)", flush=True)
        if ok and (best is None or k > best["k"]):
            n = N - k; best = {"p": p, "m": m, "r": r, "delta": delta, "k": k, "n": n, "d_lower_bound": dmin, "gamma_upper": float(np.log(n / k) / np.log(dmin)),
                               "cubic_coeffs": [int(x) for x in coef], "monomials": mons, "level": lev, "puncture_columns_1indexed": [int(x) + 1 for x in O]}
            json.dump({"best": best, "tested": tested, "cap": cap}, open(out, "w"))
    if best is not None and best["k"] >= cap - p: break
json.dump({"best": best, "tested": tested, "cap": cap}, open(out, "w"))
print("BEST:", None if best is None else {kk: v for kk, v in best.items() if kk not in ("cubic_coeffs", "monomials", "puncture_columns_1indexed")}, flush=True)
