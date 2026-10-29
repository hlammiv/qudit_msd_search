"""Weight-7 census on the class block of [[1756,260,7]]_7 (cost claim): per-output count of the
weight-7 words of ker G0 supported on class columns, by the 3+4 kernel count_weight7_per_output.
Words with non-class support are NOT counted here (see CERTIFY_P7.md); the block count is a lower
bound on L_7 and the binomial-tail bound is the rigorous upper bound.
usage: census7_1756.py OUT.json [threads-info in env]"""
import sys, json, time, numpy as np
from qmsd.weighted import weighted_support, quadric, class_points, weighted_code
from qmsd.structured_ad import _rank
rec = json.load(open("runs/weighted/p7_m4_sphere_1756_260_cols.json")); p, m, r = 7, 4, 7
GX, X, mu = weighted_support(p, m, quadric(p, m, (1, 1, 1, 1)), r)
cols = rec["cols_in_X"]; cset = set(class_points(mu, p, 1))
G0, H1, coeffs = weighted_code(GX, mu, cols, p)
keep = [l for l in range(GX.shape[1]) if l not in set(cols)]
C = [j for j, l in enumerate(keep) if l in cset]
Gc = np.ascontiguousarray(G0[:, C]); rows = [t for t in range(Gc.shape[0]) if Gc[t].any()]; Gc = np.ascontiguousarray(Gc[rows])
H1c = np.ascontiguousarray(np.asarray(H1, dtype=np.int64)[:, C] % p)
print(f"G0 {G0.shape} H1 {H1.shape}; class block {Gc.shape} rank {_rank(Gc, p)}", flush=True)
from qmsd.mindist_nb import count_weight7_per_output
t0 = time.time(); prog = []
def log(lo, hi, s, found):
    prog.append([lo, hi, round(s), found]); print(f"  cols [{lo},{hi}) {s:.0f}s, hits so far {found} (elapsed {time.time()-t0:.0f}s)", flush=True)
    json.dump({"status": "running", "progress": prog}, open(sys.argv[1], "w"))
L7, n_classes, per_output = count_weight7_per_output(Gc, H1c, p, block=8, log=log)
out = {"status": "done", "k": 260, "n": int(G0.shape[1]), "block_shape": list(Gc.shape), "weight7_block_classes": int(n_classes),
       "A7_block": int(n_classes) * (p - 1), "L7_block_per_output_mean": L7, "per_output": [int(x) for x in per_output], "seconds": round(time.time() - t0), "progress": prog}
json.dump(out, open(sys.argv[1], "w"), indent=1)
print("RESULT:", {k: v for k, v in out.items() if k not in ("per_output", "progress")}, flush=True)
