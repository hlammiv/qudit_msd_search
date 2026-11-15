"""d''(Y) >= 8 test for a paired-cap draw (probe_paired_cap.py): H = RM(6,3) shortened at Y (343 - t columns);
weights <= 5 (serial + parallel 2+3), 6 (3+3 witness) and 7 (3+4 witness).  With d'(Y) = 4 this makes the
[[686 - 2t, 2t]] single-species code d = min(2 d', d'') = 8.
usage: PYTHONPATH=. python scratchpad/cert_d2_paired.py SWEEP.json T SAMPLE OUT.json"""
import sys, json, time, numpy as np
from qmsd.reedmuller import rm_generator
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import weight_d_exists_fast, weight5_exists_parallel, weight6_witness_3x3, weight7_witness_3x4
src, t, s, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
p = 7
dr = [d for d in json.load(open(src))["draws"] if d["t"] == t and d["sample"] == s][0]
R6 = np.asarray(rm_generator(6, 3, p), dtype=np.int64) % p
H = np.asarray(shorten_matrix(R6, dr["Y"], p), dtype=np.int64) % p
res = {"t": t, "sample": s, "Y": dr["Y"], "d1": dr["d1"], "cols": int(H.shape[1])}
T0 = time.time(); d2 = None
for w in (2, 3, 4):
    if weight_d_exists_fast(H, p, w): d2 = w; break
if d2 is None and weight5_exists_parallel(H, p): d2 = 5
if d2 is None:
    w6 = weight6_witness_3x3(H, p); res["w6_witness"] = w6
    if w6 is not None: d2 = 6
print(f"weights <= 6 done: d''={d2} [{time.time()-T0:.0f}s]", flush=True)
if d2 is None:
    w7 = weight7_witness_3x4(H, p, log=lambda a, b, sec: print(f"  w7 lead cols {a}-{b} {sec:.0f}s", flush=True))
    res["w7_witness"] = w7
    if w7 is not None: d2 = 7
res["d2"] = d2 if d2 is not None else ">=8"; res["seconds"] = time.time() - T0
res["d_code"] = min(2 * dr["d1"], d2) if d2 is not None else min(2 * dr["d1"], 8)
print(res if "Y" not in res else {k: v for k, v in res.items() if k != "Y"}, flush=True)
json.dump(res, open(out, "w"), indent=1)
