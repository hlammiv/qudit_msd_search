"""Weight-7 census on the class block of a certified sphere-charge draw (p = 7, m = 4, r = 7):
per-output count of the weight-7 words of ker G0 supported on class columns, by the 3+4 kernel
count_weight7_per_output.  Complete (all weight-7 words) once scratchpad/mixed_w7_check.py has
reported no weight-7 word touching a non-class column for the same draw.
usage: census7_draw.py K SEED OUT.json      (draw from runs/weighted/plane_sweep/plane_sweep_k*.jsonl)"""
import sys, json, time, glob, numpy as np
from qmsd.weighted import weighted_support, quadric, class_points, weighted_code
from qmsd.structured_ad import _rank
k, seed, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]; p, m, r = 7, 4, 7
rec = None
for f in glob.glob("runs/weighted/plane_sweep/plane_sweep_k*.jsonl"):
    for l in open(f):
        if l.strip():
            rr = json.loads(l)
            if rr["k"] == k and rr["seed"] == seed and "cols_in_X" in rr: rec = rr
assert rec is not None
GX, X, mu = weighted_support(p, m, quadric(p, m, (1, 1, 1, 1)), r)
cols = rec["cols_in_X"]; cset = set(class_points(mu, p, 1))
G0, H1, coeffs = weighted_code(GX, mu, cols, p)
keep = [l for l in range(GX.shape[1]) if l not in set(cols)]
C = [j for j, l in enumerate(keep) if l in cset]
Gc = np.ascontiguousarray(G0[:, C]); rows = [t for t in range(Gc.shape[0]) if Gc[t].any()]; Gc = np.ascontiguousarray(Gc[rows])
H1c = np.ascontiguousarray(np.asarray(H1, dtype=np.int64)[:, C] % p)
print(f"k={k} seed={seed}: G0 {G0.shape} H1 {H1.shape}; class block {Gc.shape} rank {_rank(Gc, p)}", flush=True)
from qmsd.mindist_nb import count_weight7_per_output
t0 = time.time(); prog = []
def log(lo, hi, s, found):
    prog.append([lo, hi, round(s), found]); print(f"  cols [{lo},{hi}) {s:.0f}s, hits so far {found} (elapsed {time.time()-t0:.0f}s)", flush=True)
    json.dump({"status": "running", "k": k, "seed": seed, "progress": prog}, open(out, "w"))
L7, n_classes, per_output = count_weight7_per_output(Gc, H1c, p, block=8, log=log)
res = {"status": "done", "k": k, "seed": seed, "n": int(G0.shape[1]), "block_shape": list(Gc.shape), "weight7_block_classes": int(n_classes),
       "A7_block": int(n_classes) * (p - 1), "L7_block_per_output_mean": L7, "per_output": [int(x) for x in per_output], "seconds": round(time.time() - t0), "progress": prog}
json.dump(res, open(out, "w"), indent=1)
print("RESULT:", {kk: v for kk, v in res.items() if kk not in ("per_output", "progress")}, flush=True)
