"""Independent verification of a levelset_build.py code: (1) O is exactly a level set of the stored
polynomial (recomputed from coefficients), (2) rank of RM(r,4) on O = |O| using qmsd.structured_ad._rank
(a different elimination routine from the builder's), (3) the distance formula.
usage: levelset_verify.py CODE.json"""
import sys, json, itertools, time, numpy as np
from qmsd.reedmuller import points
from qmsd.structured_ad import _rank
b = json.load(open(sys.argv[1]))["best"]; p, m, r, delta = b["p"], b["m"], b["r"], b["delta"]; t0 = time.time()
P = np.asarray(points(m, p), dtype=np.int64); O = np.array(sorted(c - 1 for c in b["puncture_columns_1indexed"]))
POW = np.array([[pow(x, k, p) for k in range(p)] for x in range(p)], dtype=np.int64)   # exact x^k mod p
vals = np.zeros(len(P), dtype=np.int64)
for e, c in zip(b["cubic_monomials"], b["cubic_coeffs"]):
    if c: vals = (vals + c * np.prod([POW[P[:, i], e[i]] for i in range(m)], axis=0)) % p
Z = np.nonzero(vals == b["level"])[0]
print(f"level set recomputed: {len(Z)} points; equals O: {len(Z) == len(O) and (Z == O).all()}", flush=True)
mons = [e for e in itertools.product(range(p), repeat=m) if sum(e) <= r]
E = np.ones((len(O), len(mons)), dtype=np.int64)
for j, e in enumerate(mons):
    for i in range(m):
        if e[i]: E[:, j] = (E[:, j] * POW[P[O, i], e[i]]) % p
rk = _rank(E.copy(), p); print(f"rank of RM({r},4) on O by structured_ad._rank: {rk} (|O| = {len(O)}) -> {'INDEPENDENT' if rk == len(O) else 'DEPENDENT'} ({time.time()-t0:.0f}s)", flush=True)
dual = m * (p - 1) - r + delta - 1; a, bb = divmod(dual, p - 1); d = (p - bb) * p ** (m - a - 1); n = len(P) - len(O)
print(f"dual RM({dual},4): d >= {d}; n = {n}, k = {len(O)}, gamma <= {np.log(n/len(O))/np.log(d):.4f}", flush=True)
