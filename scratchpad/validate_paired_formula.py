"""Full-code check of d = min(2d', d'') on a paired-cap draw with small d': no word of weight < 2d' in ker(G0) of the
[[686 - 2t, 2t]] code, and (when 2d' = 6) a weight-6 word exists.
usage: PYTHONPATH=. python scratchpad/validate_paired_formula.py SWEEP.json T SAMPLE"""
import sys, json, time, numpy as np
from qmsd.weighted import weighted_support, product_of_linear
from qmsd.reedmuller import points
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import weight_d_exists_fast, weight5_exists_parallel, weight6_exists_via_2_4
src, t, s = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); p = 7
dr = [d for d in json.load(open(src))["draws"] if d["t"] == t and d["sample"] == s][0]
g = product_of_linear(p, 4, [((1, 0, 0, 0), c) for c in range(5)]); GX, X, mu = weighted_support(p, 4, g, 6)
P4 = np.asarray(points(4, p), dtype=np.int64)[X]; idxX = {tuple(int(v) for v in x): i for i, x in enumerate(P4)}
P3 = np.asarray(points(3, p), dtype=np.int64)
S = sorted(idxX[(c,) + tuple(int(v) for v in P3[y])] for y in dr["Y"] for c in (5, 6))
G0 = np.asarray(shorten_matrix(GX, S, p), dtype=np.int64) % p
T0 = time.time(); found = next((w for w in (2, 3, 4) if weight_d_exists_fast(G0, p, w)), None)
if found is None and weight5_exists_parallel(G0, p): found = 5
w6 = bool(weight6_exists_via_2_4(G0, p)) if found is None else None
print(f"t={t} s={s} d'={dr['d1']}: predicted d = {2 * dr['d1']}; full code min weight <= 5: {found}; weight 6 exists: {w6} [{time.time()-T0:.0f}s]", flush=True)
