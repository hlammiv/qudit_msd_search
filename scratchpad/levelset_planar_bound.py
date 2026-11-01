"""Exhaustive planar upper bound for a level-set code: over ALL planes of F_p^4 and all parallel
classes of lines in each, the 5-line (p=11: (p-1)/2 lines) minimum-weight dual word with the most
output points gives the smallest planar kernel word.  usage: levelset_planar_bound.py CODE.json OUT.json"""
import sys, json, itertools, time, numpy as np
from qmsd.reedmuller import points
src = json.load(open(sys.argv[1])); rec = src.get("best", src)
p, m = rec["p"], rec["m"]; O = sorted(c - 1 for c in rec["puncture_columns_1indexed"]); t0 = time.time()
P = np.asarray(points(m, p), dtype=np.int64); N = len(P); inO = np.zeros(N, dtype=bool); inO[O] = True
# dual degree 2(p-1)+b: min-weight words are (p-b) parallel lines in a plane
dual = m * (p - 1) - rec["r"] - 1; a, b = divmod(dual, p - 1); nl = p - b
# planes = fibres of surjective linear maps F^4 -> F^2 (rref 2x4), lines inside = fibres of a further map to F
normals = [v for v in itertools.product(range(p), repeat=m) if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1]
def rref2(a_, b_):
    M = np.array([a_, b_]) % p
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
print(f"{len(subs)} plane-quotient maps; dual RM({dual},4): min-weight words = {nl} parallel lines ({time.time()-t0:.0f}s)", flush=True)
best = (0, None)
# a plane-quotient map S (2x4) has fibres = planes; inside a plane the lines parallel to a direction u in ker S are the
# fibres of a third functional f with (S, f) of rank 3.  Enumerate f over normals with rank 3 and group by (S-value, f-value).
Ocnt_cache = {}
for S in subs:
    key2 = (P @ S.T) % p; key2 = key2[:, 0] * p + key2[:, 1]
    for f in normals:
        if np.linalg.matrix_rank(np.vstack([S, f]).astype(float)) < 3: continue   # fast approximate rank test, then exact below
        M3 = np.vstack([S, f]) % p
        # exact rank 3 check mod p
        from qmsd.structured_ad import _rank
        if _rank(M3.copy(), p) < 3: continue
        key3 = key2 * p + (P @ np.array(f)) % p
        cnt = np.bincount(key3[inO], minlength=p ** 3).reshape(p * p, p)   # [plane, line] -> outputs on that line
        top = -np.sort(-cnt, axis=1)[:, :nl].sum(axis=1); mx = int(top.max())
        if mx > best[0]: best = (mx, (S.tolist(), list(f), int(np.argmax(top)))); print(f"  new max outputs on {nl} parallel lines: {mx} -> kernel weight {nl*p - mx} ({time.time()-t0:.0f}s)", flush=True)
json.dump({"max_outputs_on_minweight_support": best[0], "kernel_weight_upper": nl * p - best[0], "where": best[1], "seconds": round(time.time() - t0)}, open(sys.argv[2], "w"))
print("RESULT: d <=", nl * p - best[0], f"({time.time()-t0:.0f}s)", flush=True)
