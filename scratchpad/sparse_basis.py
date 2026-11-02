"""Lighter Clifford circuits: (a) a sparse basis for the checks (rowspace of G0) from structured light
vectors + greedy selection, completed by row-reduced rows; (b) light representatives of each output's
logical operator (H1_i + rowspace G0), by greedy descent: repeatedly add the multiple of a check
(light basis or row-reduced) that lowers the weight most.  Reports gates for checks and decoding
before/after under the same model as clifford_cost.py.
usage: sparse_basis.py ARTIFACT_ID OUT.json [n_logical_samples]"""
import sys, json, time, itertools, numpy as np, numba as nb
from qmsd.results import load_result_catalog
from qmsd.performance import record_matrices
from qmsd.reedmuller import points
aid, out = sys.argv[1], sys.argv[2]; nsamp = int(sys.argv[3]) if len(sys.argv) > 3 else 40
r = next(x for x in load_result_catalog() if x.artifact_id == aid); p = r.p; t0 = time.time()
G0, H1 = (__import__('scratchpad.fast_matrices', fromlist=['x']).fast_matrices(r) if r.n > 5000 and (r.construction or {}).get('family') == 'punctured_rm' else record_matrices(r)); G0 = np.asarray(G0, dtype=np.int64) % p; H1 = np.asarray(H1, dtype=np.int64) % p
G0 = G0[[i for i in range(G0.shape[0]) if G0[i].any()]]; n = G0.shape[1]; k = H1.shape[0]

@nb.njit(cache=True)
def rref_mod(A, p):
    A = A.copy(); R, C = A.shape; inv = np.zeros(p, dtype=np.int64)
    for x in range(1, p):
        for y in range(1, p):
            if (x * y) % p == 1: inv[x] = y
    piv = np.full(R, -1, dtype=np.int64); rk = 0
    for c in range(C):
        pr = -1
        for i in range(rk, R):
            if A[i, c] != 0: pr = i; break
        if pr < 0: continue
        if pr != rk:
            for j in range(C): t = A[pr, j]; A[pr, j] = A[rk, j]; A[rk, j] = t
        iv = inv[A[rk, c]]
        for j in range(C): A[rk, j] = (A[rk, j] * iv) % p
        for i in range(R):
            if i != rk and A[i, c] != 0:
                f = A[i, c]
                for j in range(C): A[i, j] = (A[i, j] - f * A[rk, j]) % p
        piv[rk] = c; rk += 1
        if rk == R: break
    return A[:rk], piv[:rk]

@nb.njit(cache=True)
def greedy_descend(v, B, p, passes):
    """Lower wt(v) by adding multiples of rows of B while it helps."""
    v = v.copy(); R, C = B.shape
    w = 0
    for j in range(C):
        if v[j] != 0: w += 1
    for _ in range(passes):
        improved = False
        for i in range(R):
            bestc = 0; bestw = w
            for c in range(1, p):
                ww = 0
                for j in range(C):
                    if (v[j] + c * B[i, j]) % p != 0: ww += 1
                if ww < bestw: bestw = ww; bestc = c
            if bestc:
                for j in range(C): v[j] = (v[j] + bestc * B[i, j]) % p
                w = bestw; improved = True
        if not improved: break
    return v, w

Gs, piv = rref_mod(G0, p); rk = Gs.shape[0]
w_rref = (Gs != 0).sum(axis=1)
print(f"{r.label}: G0 rank {rk}, rref check weights mean {w_rref.mean():.0f} ({time.time()-t0:.0f}s)", flush=True)

# ---- (a) structured light check vectors, code-family specific
cands = []
fam = (r.construction or {}).get("family")
if "level" in (r.note or "") and fam == "punctured_rm":
    # g * (hyperplane indicator): g is the cubic whose zero set is O; nonzero off O, so support = H minus O
    meta = json.load(open("scratchpad/routes/algebraic/p11_cubic_best.json")) if p == 11 else json.load(open(f"runs/rm_m4/p{p}_levelset.json"))["best"]
    P = np.asarray(points(4, p), dtype=np.int64); O = sorted(c - 1 for c in r.puncture_columns); keep = np.array([j for j in range(len(P)) if j not in set(O)])
    mons = [e for e in itertools.product(range(4), repeat=4) if sum(e) <= 3]
    coef = meta.get("coeffs") or meta.get("cubic_coeffs")
    gv = np.zeros(len(P), dtype=np.int64)
    for e, c in zip(mons, coef):
        if c: gv = (gv + c * np.prod([np.power(P[:, i], e[i]) % p for i in range(4)], axis=0)) % p
    gv = (gv - meta["level"]) % p
    normals = [v for v in itertools.product(range(p), repeat=4) if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1]
    rng = np.random.default_rng(0)
    for v in normals:
        for c0 in range(p):
            ind = (((P @ np.array(v)) % p) == c0).astype(np.int64)
            cands.append(((gv * ind) % p)[keep])
            if len(cands) >= 6 * rk: break
        if len(cands) >= 6 * rk: break
elif fam == "weighted_rm":
    # (g-1)(g-6) * (product of three linear forms): vanishes on the class block and on three hyperplanes
    from qmsd.weighted import weighted_support, quadric
    GX, X, mu = weighted_support(7, 4, quadric(7, 4, (1, 1, 1, 1)), 7)
    P = np.asarray(points(4, 7), dtype=np.int64); outs = set(r.construction["outputs_in_X"])
    keepX = [l for l in range(len(X)) if l not in outs]; Pk = P[np.array([int(X[l]) for l in keepX])]
    gq = (Pk ** 2).sum(axis=1) % 7; base = ((gq - 1) * (gq - 6)) % 7
    # the stored code rescales columns by charges; recover the column scaling from G0 itself by projection below
    normals = [v for v in itertools.product(range(7), repeat=4) if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1]
    rng = np.random.default_rng(0)
    for _ in range(8 * rk):
        f = base.copy()
        for _j in range(3):
            v = np.array(normals[rng.integers(len(normals))]); c0 = rng.integers(7)
            f = (f * ((Pk @ v - c0) % 7)) % 7
        cands.append(f)
cands = [c for c in cands if c.any()]
print(f"{len(cands)} structured candidate vectors ({time.time()-t0:.0f}s)", flush=True)
# keep only candidates in rowspace(G0) (the weighted code's columns carry charge scalings; test membership)
def in_rowspace(v):
    vv = v.copy() % p
    for i, c in enumerate(piv):
        if vv[c]: vv = (vv - vv[c] * Gs[i]) % p
    return not vv.any()
ech = {}                           # incremental echelon form: pivot column -> normalised row
def add_if_independent(v):
    vv = v.copy() % p
    for c, row in ech.items():
        if vv[c]: vv = (vv - vv[c] * row) % p
    nz = np.nonzero(vv)[0]
    if len(nz) == 0: return False
    c = int(nz[0]); vv = (vv * pow(int(vv[c]), -1, p)) % p
    for c2 in list(ech):
        if ech[c2][c]: ech[c2] = (ech[c2] - ech[c2][c] * vv) % p
    ech[c] = vv; return True
chosen = []
if cands:
    ok = [c for c in cands[: 2 * rk] if in_rowspace(c)]
    print(f"  {len(ok)} of {min(len(cands), 2*rk)} candidates lie in the check space", flush=True)
    ok.sort(key=lambda c: int((c != 0).sum()))
    for c in ok:
        if add_if_independent(c): chosen.append(c)
        if len(chosen) == rk: break
nl = len(chosen)
for row in Gs:
    if len(chosen) == rk: break
    if add_if_independent(row): chosen.append(row)
light_full = np.array(chosen)
w_light = (light_full != 0).sum(axis=1)
print(f"sparse check basis: {nl} structured + {rk - nl} completed; check gates {int(w_rref.sum())} -> {int(w_light.sum())} ({w_rref.sum()/w_light.sum():.2f}x) ({time.time()-t0:.0f}s)", flush=True)

# ---- (b) light logical representatives (sampled outputs)
def reduce_vs(H, Gs, piv):
    H = H.copy() % p
    for i, c in enumerate(piv): H = (H - np.outer(H[:, c], Gs[i])) % p
    return H
rng = np.random.default_rng(1); samp = rng.choice(k, size=min(nsamp, k), replace=False)
Hr = reduce_vs(H1[samp], Gs, piv)
w_before = (Hr != 0).sum(axis=1)
B = np.ascontiguousarray(np.vstack([light_full, Gs]))
w_after = []
for i in range(len(samp)):
    v, w = greedy_descend(np.ascontiguousarray(Hr[i]), B, p, 3); w_after.append(w)
w_after = np.array(w_after)
print(f"logical reps ({len(samp)} sampled outputs): mean weight {w_before.mean():.0f} -> {w_after.mean():.0f} ({w_before.mean()/w_after.mean():.2f}x), min {w_after.min()} ({time.time()-t0:.0f}s)", flush=True)
dec_before = w_before.mean() * k; dec_after = w_after.mean() * k
tot_b = w_rref.sum() + dec_before; tot_a = w_light.sum() + dec_after
res = {"code": r.label, "n": n, "k": k, "checks": rk, "structured_light_rows": nl, "check_gates_rref": int(w_rref.sum()), "check_gates_sparse": int(w_light.sum()),
       "logical_weight_rref_mean": float(w_before.mean()), "logical_weight_descended_mean": float(w_after.mean()), "logical_weight_descended_min": int(w_after.min()),
       "decode_gates_rref_est": float(dec_before), "decode_gates_sparse_est": float(dec_after), "gates_per_output_before": float(tot_b / k), "gates_per_output_after": float(tot_a / k),
       "saving": float(tot_b / tot_a), "seconds": round(time.time() - t0)}
json.dump(res, open(out, "w"), indent=1); print(json.dumps(res), flush=True)
