"""Explicit low-weight witnesses for sweep survivors, from flats only (cheap, fits a 2 GB box).
For each survivor (d_plane >= 7, no block word <= 5): re-scan the affine planes and record the
lightest planar word (G0 columns, coefficients, points), verified against the full G0. If
d_plane >= 8, scan the 2800 affine 3-flats with the 3+4 weight-7 witness kernel. A witness
proves d <= weight; the lower bound still needs the full scans.
usage: flat_witness.py OUT.jsonl K [SEED ...]   (no seeds = all survivors of that k)"""
import sys, json, glob, itertools, time, os, numpy as np
from qmsd.reedmuller import points
from qmsd.weighted import weighted_support, quadric, class_points
from qmsd.structured_ad import _rank, _null_space
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import weight_d_exists_fast, weight7_witness_3x4
out, k = sys.argv[1], int(sys.argv[2]); seeds = [int(s) for s in sys.argv[3:]]
p, m, r = 7, 4, 7
P = np.asarray(points(m, p), dtype=np.int64)
g = quadric(p, m, (1, 1, 1, 1)); GX, X, mu = weighted_support(p, m, g, r); cand = class_points(mu, p, 1)
normals = np.array([v for v in itertools.product(range(p), repeat=m) if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1])
def rref2(a, b):
    M = np.array([a, b]) % p
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
Xc = np.array([int(X[l]) for l in cand]); PX = P[Xc]
COEFS = {nu: np.array(list(itertools.product(range(p), repeat=nu)))[1:] for nu in range(1, 7)}
def groups(keyfun_rows):
    """yield index arrays of class points sharing a key (one flat each)."""
    key = keyfun_rows; order = np.argsort(key, kind="stable"); ks = key[order]
    starts = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]]); ends = np.r_[starts[1:], len(ks)]
    for a, b in zip(starts, ends): yield order[a:b]
recs = []
for f in glob.glob(f"runs/weighted/plane_sweep/plane_sweep_k{k}*.jsonl"):
    for l in open(f):
        if l.strip():
            rr = json.loads(l)
            if rr["d_plane"] >= 7 and rr.get("block_d_le5") is None and (not seeds or rr["seed"] in seeds): recs.append(rr)
done = set()
if os.path.exists(out):
    for l in open(out):
        if l.strip(): done.add(json.loads(l)["seed"])
for rr in sorted(recs, key=lambda x: x["seed"]):
    if rr["seed"] in done: continue
    cols = rr["cols_in_X"]; t0 = time.time()
    G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p
    keep = [l for l in range(GX.shape[1]) if l not in set(cols)]; pos = {l: j for j, l in enumerate(keep)}
    sel = set(cols); surv = np.array([cand[i] not in sel for i in range(len(cand))])
    colidx = np.array([pos.get(cand[i], -1) for i in range(len(cand))])
    local = -np.ones(len(cand), dtype=int); local[np.nonzero(surv)[0]] = np.arange(surv.sum())
    cls_cols = colidx[surv]; Gc = G0[:, cls_cols]; Gc = Gc[[t for t in range(Gc.shape[0]) if Gc[t].any()]]
    best = (99, None)
    for Sb in subs:
        v = (PX @ Sb.T) % p
        for idx in groups(v[:, 0] * p + v[:, 1]):
            cc = local[idx]; cc = cc[cc >= 0]
            if len(cc) < 2: continue
            M = Gc[:, cc]
            if _rank(M.T.copy(), p) == len(cc): continue
            ns = np.asarray(_null_space(M, p)) % p; nu = ns.shape[0]
            if nu > 6: continue
            W = (COEFS[nu] @ ns) % p; wts = (W != 0).sum(axis=1); j = int(wts.argmin())
            if wts[j] < best[0]: best = (int(wts[j]), (cc, W[j]))
    rec = {"k": k, "seed": rr["seed"], "n": G0.shape[1], "d_plane": best[0]}
    if best[1] is not None:
        cc, w = best[1]; nz = np.flatnonzero(w)
        g0cols = [int(cls_cols[c]) for c in cc[nz]]; coef = [int(w[i]) for i in nz]
        vec = np.zeros(G0.shape[1], dtype=np.int64); vec[g0cols] = coef
        assert ((G0 @ vec) % p == 0).all()
        rec.update({"witness_weight": len(g0cols), "witness_G0_cols": g0cols, "witness_coeffs": coef,
                    "witness_points": [tuple(int(x) for x in P[int(X[keep[c]])]) for c in g0cols], "witness_from": "plane"})
    if best[0] >= 8:
        # 3-flats: cosets of the 400 hyperplane directions ... a 3-flat is {x : a.x = c}
        found = None; t1 = time.time(); n3 = 0
        for a in normals:
            for idx in groups((PX @ a) % p):
                cc = local[idx]; cc = cc[cc >= 0]
                if len(cc) < 7: continue
                n3 += 1; M = np.ascontiguousarray(Gc[:, cc])
                if any(weight_d_exists_fast(M, p, d) for d in (2, 3, 4, 5, 6)):
                    found = ("<=6", cc, None); break
                w7 = weight7_witness_3x4(M, p)
                if w7 is not None: found = ("7", cc, w7); break
            if found: break
        rec["flats3_scanned"] = n3; rec["t_flats3"] = round(time.time() - t1)
        if found is None:
            rec["d_flat3"] = ">=8"
        else:
            tag, cc, w7 = found
            if w7 is not None:
                g0cols = [int(cls_cols[cc[c]]) for c in w7]
                ns = np.asarray(_null_space(G0[:, g0cols], p)) % p; assert ns.shape[0] >= 1
                coef = [int(x) for x in ns[0]]
                vec = np.zeros(G0.shape[1], dtype=np.int64); vec[g0cols] = coef; assert ((G0 @ vec) % p == 0).all()
                rec.update({"d_flat3": 7, "witness_weight": 7, "witness_G0_cols": g0cols, "witness_coeffs": coef,
                            "witness_points": [tuple(int(x) for x in P[int(X[keep[c]])]) for c in g0cols], "witness_from": "3-flat"})
            else:
                rec["d_flat3"] = tag
    rec["t_total"] = round(time.time() - t0)
    with open(out, "a") as f: f.write(json.dumps(rec) + "\n")
    print(json.dumps({a: b for a, b in rec.items() if a not in ("witness_G0_cols",)}), flush=True)
