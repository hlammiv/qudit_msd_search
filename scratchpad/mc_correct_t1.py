"""Monte Carlo of a detect-and-correct round with the t = 1 lookup decoder on a catalog code:
syndrome 0 -> keep; syndrome = a * (column j of G0) for a unique (j, a) -> undo Z^a on qudit j, keep;
else discard.  Reports acceptance and per-output / block logical error for t = 0 and t = 1 on the same shots.
usage: mc_correct_t1.py ARTIFACT_ID EPS SHOTS OUT.json"""
import sys, json, time, numpy as np
from qmsd.results import load_result_catalog
from qmsd.performance import record_matrices
aid, eps, shots, out = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
r = next(x for x in load_result_catalog() if x.artifact_id == aid); p = r.p
G0, H1 = record_matrices(r); G0 = np.asarray(G0, dtype=np.int64) % p; H1 = np.asarray(H1, dtype=np.int64) % p
G0 = G0[[t for t in range(G0.shape[0]) if G0[t].any()]]; m, n = G0.shape; k = H1.shape[0]
pw = np.array([pow(p, i, 2 ** 62) for i in range(m)], dtype=np.int64)  # hash of a syndrome vector (collision-checked below)
def key(S): return (S * pw).sum(axis=-1)
table = {}
amb = 0
for j in range(n):
    for a in range(1, p):
        kk = int(key((a * G0[:, j]) % p))
        if kk in table: amb += 1; table[kk] = None
        else: table[kk] = (j, a)
print(f"{r.label}: G0 {G0.shape}, lookup table {len(table)} keys, {amb} ambiguous (dropped)", flush=True)
rng = np.random.default_rng(7); t0 = time.time(); B = 10_000
stats = {"shots": 0, "acc0": 0, "acc1": 0, "fail_out0": 0, "fail_blk0": 0, "fail_out1": 0, "fail_blk1": 0, "corrected": 0}
while stats["shots"] < shots:
    b = min(B, shots - stats["shots"])
    e = np.where(rng.random((b, n)) < eps, rng.integers(1, p, size=(b, n)), 0)
    S = (e @ G0.T) % p; z = ~S.any(axis=1)
    L = (e @ H1.T) % p
    stats["acc0"] += int(z.sum()); f0 = L[z] != 0; stats["fail_out0"] += int(f0.sum()); stats["fail_blk0"] += int(f0.any(axis=1).sum())
    nz = np.nonzero(~z)[0]; ks = key(S[nz])
    kept = 0; fo = 0; fb = 0
    for idx, kv in zip(nz, ks):
        hit = table.get(int(kv))
        if not hit: continue
        j, a = hit
        if not (((a * G0[:, j]) % p) == S[idx]).all(): continue   # hash collision guard
        ec = e[idx].copy(); ec[j] = (ec[j] - a) % p
        Lc = (ec @ H1.T) % p; kept += 1; fo += int((Lc != 0).sum()); fb += int((Lc != 0).any())
    stats["acc1"] += int(z.sum()) + kept; stats["corrected"] += kept
    stats["fail_out1"] += int(f0.sum()) + fo; stats["fail_blk1"] += int(f0.any(axis=1).sum()) + fb
    stats["shots"] += b
s = stats; res = {"code": r.label, "eps": eps, **s, "P_acc_t0": s["acc0"] / s["shots"], "P_acc_t1": s["acc1"] / s["shots"], "P_acc_t0_theory": (1 - eps) ** n,
                  "P_acc_t1_lower_theory": (1 - eps) ** n + n * eps * (1 - eps) ** (n - 1),
                  "eps_out_t0": s["fail_out0"] / max(1, k * s["acc0"]), "eps_out_t1": s["fail_out1"] / max(1, k * s["acc1"]),
                  "block_t0": s["fail_blk0"] / max(1, s["acc0"]), "block_t1": s["fail_blk1"] / max(1, s["acc1"]), "seconds": round(time.time() - t0)}
json.dump(res, open(out, "w"), indent=1); print(json.dumps(res), flush=True)
