"""Geometry of the low-weight words found on the sphere-charge codes: affine span, collinear
triples, coplanar quadruples, which sphere (g = 1 or 6) each point is on, and for the d = 2 / 4
sets (k = 290 / 280, seed 2026) and the d = 3 set (k = 265 sample 0, seed 7) the actual words
from the streamed MITM."""
import json, random, itertools, numpy as np
from qmsd.reedmuller import points
from qmsd.weighted import weighted_support, quadric, class_points
from qmsd.structured_ad import _rank, mitm_streamed
from qmsd.triorthogonal import shorten_matrix
p, m, r = 7, 4, 7
P = np.asarray(points(m, p), dtype=np.int64)
g = quadric(p, m, (1, 1, 1, 1)); GX, X, mu = weighted_support(p, m, g, r); cand = class_points(mu, p, 1)
def sample(seed, klist):
    rng = random.Random(seed); out = {}
    for k in klist:
        cols = sorted(rng.sample(cand, k))
        while _rank(GX[:, cols].T.copy(), p) < k: cols = sorted(rng.sample(cand, k))
        out[k] = cols
    return out
def describe(pts, vals=None):
    Q = P[pts]; A = Q[1:] - Q[0]
    span = _rank(A % p, p) if len(pts) > 1 else 0
    col = sum(1 for t in itertools.combinations(range(len(pts)), 3) if _rank((Q[list(t)][1:] - Q[t[0]]) % p, p) <= 1)
    cop = sum(1 for t in itertools.combinations(range(len(pts)), 4) if _rank((Q[list(t)][1:] - Q[t[0]]) % p, p) <= 2)
    sph = [int(g[x]) for x in pts]
    return f"span {span}, collinear triples {col}, coplanar quads {cop}, g-values {sph}" + (f", coeffs {vals}" if vals is not None else "")
def words(cols, d_max):
    G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p
    keep = [l for l in range(GX.shape[1]) if l not in set(cols)]
    for d in range(2, d_max + 1):
        ws = list(itertools.islice(mitm_streamed(G0, d, p), 6))
        if ws:
            print(f"   weight-{d} words (first {len(ws)} scalar classes):")
            for supp, vals in ws:
                pts = [int(X[keep[c]]) for c in supp]
                print("     ", [tuple(int(v) for v in P[x]) for x in pts], "|", describe(pts, list(map(int, vals))))
            return
# the certified witnesses
rec270 = json.load(open("runs/weighted/p7_m4_sphere_1746_270_cols.json")); rec260 = json.load(open("runs/weighted/p7_m4_sphere_1756_260_cols.json"))
for name, rec, key, cf in (("k=270 weight-6", rec270, "weight6_witness_G0_cols", "witness_coeffs"), ("k=260 weight-7", rec260, "weight7_witness_G0_cols", "witness7_coeffs")):
    keep = [l for l in range(GX.shape[1]) if l not in set(rec["cols_in_X"])]
    pts = [int(X[keep[c]]) for c in rec[key]]
    print(name, [tuple(int(v) for v in P[x]) for x in pts]); print("   ", describe(pts, rec[cf][0]))
# the weak draws
s2026 = sample(2026, [290, 280]); s7 = sample(7, [275, 275, 265])   # seed 7: k=275 x2 then k=265 sample 0
print("k=290 (d=2):"); words(s2026[290], 2)
print("k=280 (d=4):"); words(s2026[280], 4)
