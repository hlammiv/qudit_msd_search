"""Distance stages for an rm_m4_candidate.py output: weight <= 4 by the serial numba MITM, weight 5 by the
thread-parallel 2+3 scan.  Resumable; appends results into the candidate json.
usage: rm_m4_scan.py CAND.json [--w5]"""
import sys, json, time, numpy as np
from qmsd.reedmuller import rm_generator
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import weight_d_exists_fast, weight5_exists_parallel
f = sys.argv[1]; W5 = "--w5" in sys.argv
c = json.load(open(f)); p, m, r = c["p"], c["m"], c["r"]; cols0 = [x - 1 for x in c["puncture_columns_1indexed"]]
G = np.asarray(rm_generator(r, m, p), dtype=np.int64) % p
G0 = np.ascontiguousarray(np.asarray(shorten_matrix(G, cols0, p), dtype=np.int64) % p)
G0 = G0[[t for t in range(G0.shape[0]) if G0[t].any()]]
print(f"G0 {G0.shape} (dim A = {c['dim_A']})", flush=True)
c.setdefault("stages", {})
for d in (2, 3, 4):
    if str(d) in c["stages"]: continue
    t = time.time(); ex = bool(weight_d_exists_fast(G0, p, d)); c["stages"][str(d)] = {"exists": ex, "seconds": round(time.time() - t)}
    json.dump(c, open(f, "w"), indent=1); print(f"weight {d} exists: {ex} ({time.time()-t:.0f}s)", flush=True)
    if ex: c["d"] = d; json.dump(c, open(f, "w"), indent=1); sys.exit(0)
c["d_lower"] = 5; json.dump(c, open(f, "w"), indent=1); print("d >= 5", flush=True)
if W5 and "5" not in c["stages"]:
    t = time.time(); ex = bool(weight5_exists_parallel(G0, p)); c["stages"]["5"] = {"exists": ex, "seconds": round(time.time() - t)}
    c["d" if ex else "d_lower"] = 5 if ex else 6; json.dump(c, open(f, "w"), indent=1); print(f"weight 5 exists: {ex} ({time.time()-t:.0f}s)", flush=True)
