"""Weight-6 test for draws of probe_2hp.py that cleared weight 5: rebuild G0 from the stored output columns and
run the 3+3 witness (left table ~20 GB at n ~ 570: lenore only).  Writes OUT.json with witness or None.
usage: PYTHONPATH=. python scratchpad/cert6_2hp.py PROBE.json K SAMPLE OUT.json"""
import sys, json, time, numpy as np
from qmsd.weighted import weighted_support, product_of_linear
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import weight_d_exists_fast, weight5_exists_parallel, weight6_witness_3x3
src, k, s, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
p, m, r = 7, 4, 6
g = product_of_linear(p, m, [((1, 0, 0, 0), c) for c in range(5)])
GX, X, mu = weighted_support(p, m, g, r)
dr = [d for d in json.load(open(src))["draws"] if d["k"] == k and d["sample"] == s][0]
G0 = np.asarray(shorten_matrix(GX, dr["cols_in_X"], p), dtype=np.int64) % p
assert not any(weight_d_exists_fast(G0, p, w) for w in (2, 3, 4)) and not weight5_exists_parallel(G0, p)
t0 = time.time(); w = weight6_witness_3x3(G0, p)
res = {"k": k, "sample": s, "n": int(G0.shape[1]), "weight6_witness": None if w is None else [int(c) for c in w],
       "d": 6 if w is not None else ">=7", "seconds": time.time() - t0}
print(res, flush=True); json.dump(res, open(out, "w"), indent=1)
