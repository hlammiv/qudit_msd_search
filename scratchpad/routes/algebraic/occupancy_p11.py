"""Occupancy statistics of the k=2290 output set of runs/rm_m4/p11_r13_k2290_seed1.json:
max outputs on any affine line / plane / hyperplane (and the direction-wise line maxima).
Run: ulimit -v 4000000; OMP_NUM_THREADS=1 PYTHONPATH=. python scratchpad/routes/algebraic/occupancy_p11.py"""
import json, sys, time
import numpy as np
from itertools import product
from qmsd.reedmuller import points

p, m = 11, 4
d = json.load(open("runs/rm_m4/p11_r13_k2290_seed1.json"))
O = np.array(d["puncture_columns_1indexed"]) - 1
pts = points(m, p)                       # (p^m, m), x_1 least significant
n0 = p ** m
isO = np.zeros(n0, bool); isO[O] = True
k = len(O)
print("k =", k, "n =", n0 - k)

def encode(X):
    return (X % p) @ (p ** np.arange(m))

# --- directions (projective points): canonical reps with first nonzero coord = 1
def directions(j):
    """All j-dim linear subspaces? For j=1 just direction vectors."""
    reps = []
    for v in product(range(p), repeat=m):
        v = np.array(v)
        nz = np.nonzero(v)[0]
        if len(nz) and v[nz[0]] == 1:
            reps.append(v)
    return np.array(reps)

t0 = time.time()
dirs = directions(1)
assert len(dirs) == (p**m - 1)//(p - 1)
line_max = 0; hist = np.zeros(p + 1, int)
for v in dirs:
    i = np.nonzero(v)[0][0]
    base = (pts - np.outer(pts[:, i], v)) % p      # representative with x_i = 0
    lab = encode(base)
    cnt = np.bincount(lab[isO], minlength=n0)
    cnt = cnt[np.bincount(lab, minlength=n0) > 0]  # only real lines
    hist += np.bincount(cnt, minlength=p + 1)
    line_max = max(line_max, cnt.max())
print(f"lines: max occupancy {line_max}; histogram (occ: #lines) =",
      {i: int(hist[i]) for i in range(p + 1) if hist[i]}, f"({time.time()-t0:.1f}s)")

# --- hyperplanes: a.x = c for a in dirs (as normals), c in F_p
hyp_max = 0; hyp_hist = {}
for a in dirs:
    val = (pts @ a) % p
    cnt = np.bincount(val[isO], minlength=p)
    hyp_max = max(hyp_max, cnt.max())
    for c in cnt: hyp_hist[int(c)] = hyp_hist.get(int(c), 0) + 1
print("hyperplanes: max occupancy", hyp_max, "mean", k * p**(m-1) / n0)

# --- planes: 2-dim subspaces U (pairs of directions, canonical), cosets labelled by
# the projection onto a complement.  Enumerate U by reduced row echelon bases.
t0 = time.time()
plane_max = 0; nplanes = 0; plane_hist = {}
seen = set()
# RREF 2x4 matrices over F_p: pivot columns (i<j)
for i in range(m):
    for j in range(i + 1, m):
        free_i = [c for c in range(i + 1, m) if c != j]
        free_j = [c for c in range(j + 1, m)]
        for vi in product(range(p), repeat=len(free_i)):
            for vj in product(range(p), repeat=len(free_j)):
                u = np.zeros(m, int); w = np.zeros(m, int)
                u[i] = 1; w[j] = 1
                for c, val in zip(free_i, vi): u[c] = val
                for c, val in zip(free_j, vj): w[c] = val
                # coset label: x - x_i u - (x_j - x_i u_j) w  has coords i and j zero
                base = (pts - np.outer(pts[:, i], u)) % p
                base = (base - np.outer(base[:, j], w)) % p
                lab = encode(base)
                tot = np.bincount(lab, minlength=n0)
                cnt = np.bincount(lab[isO], minlength=n0)[tot > 0]
                assert (tot[tot > 0] == p * p).all()
                nplanes += len(cnt)
                plane_max = max(plane_max, cnt.max())
                for c in np.unique(cnt): plane_hist[int(c)] = plane_hist.get(int(c), 0) + int((cnt == c).sum())
print(f"planes: {nplanes} planes, max occupancy {plane_max}, mean {k*p*p/n0:.2f} ({time.time()-t0:.1f}s)")
print("plane occupancy tail:", {c: plane_hist[c] for c in sorted(plane_hist) if c >= 35})
print("hyperplane occupancy tail:", {c: hyp_hist[c] for c in sorted(hyp_hist) if c >= 250})
