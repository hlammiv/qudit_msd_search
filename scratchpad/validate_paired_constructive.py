"""Constructive check of d = 2d' on a paired-cap draw with d' = 3: find a weight-3 word a of RM(5,3) shortened at Y,
lift it to (a | -a) on the two hyperplanes off S and verify it lies in ker(G0) (weight 6); then confirm no word of
weight <= 5 in the full code.  usage: PYTHONPATH=. python scratchpad/validate_paired_constructive.py SWEEP.json T SAMPLE"""
import sys, json, time, itertools, numpy as np
from qmsd.weighted import weighted_support, product_of_linear
from qmsd.reedmuller import points, rm_generator
from qmsd.triorthogonal import shorten_matrix
from qmsd.structured_ad import _null_space
from qmsd.pointsets import batched_rank_modp
from qmsd.mindist_nb import weight_d_exists_fast, weight5_exists_parallel
src, t, s = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); p = 7
dr = [d for d in json.load(open(src))["draws"] if d["t"] == t and d["sample"] == s][0]; Y = dr["Y"]
P3 = np.asarray(points(3, p), dtype=np.int64); R5 = np.asarray(rm_generator(5, 3, p), dtype=np.int64) % p
offY = [x for x in range(p ** 3) if x not in set(Y)]
H3 = np.asarray(shorten_matrix(R5, Y, p), dtype=np.int64) % p          # columns = offY in order
tr = np.array(list(itertools.combinations(range(len(offY)), 3)))
w3 = None
for i in range(0, len(tr), 200000):
    blk = tr[i:i + 200000]; rk = batched_rank_modp(np.transpose(H3[:, blk], (1, 0, 2)), p)
    hit = np.nonzero(rk < 3)[0]
    if len(hit): w3 = blk[hit[0]]; break
ns = np.asarray(_null_space(H3[:, list(w3)], p)) % p; coef = ns[0]
assert (coef != 0).all()
g = product_of_linear(p, 4, [((1, 0, 0, 0), c) for c in range(5)]); GX, X, mu = weighted_support(p, 4, g, 6)
P4 = np.asarray(points(4, p), dtype=np.int64)[X]; idxX = {tuple(int(v) for v in x): i for i, x in enumerate(P4)}
S = sorted(idxX[(c,) + tuple(int(v) for v in P3[y])] for y in Y for c in (5, 6)); keep = [l for l in range(len(X)) if l not in set(S)]
G0 = np.asarray(shorten_matrix(GX, S, p), dtype=np.int64) % p
e = np.zeros(len(keep), dtype=np.int64); pos = {l: j for j, l in enumerate(keep)}
for c3, a in zip(w3, coef):
    y = tuple(int(v) for v in P3[offY[c3]])
    e[pos[idxX[(5,) + y]]] = a % p; e[pos[idxX[(6,) + y]]] = (-a) % p
in_ker = not ((G0 @ e) % p).any()
pts = [tuple(int(v) for v in P3[offY[c]]) for c in w3]
print(f"t={t} s={s}: weight-3 d' word at y = {pts} coeffs {coef.tolist()}; lifted weight-6 word in ker G0: {in_ker}", flush=True)
T0 = time.time(); low = next((w for w in (2, 3, 4) if weight_d_exists_fast(G0, p, w)), None)
if low is None and weight5_exists_parallel(G0, p): low = 5
print(f"full code word of weight <= 5: {low}  -> d = {6 if in_ker and low is None else '?'} (predicted 2d' = 6) [{time.time()-T0:.0f}s]", flush=True)
