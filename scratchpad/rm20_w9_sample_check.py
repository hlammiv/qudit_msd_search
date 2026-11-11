"""Sample check of "every weight-8 and weight-9 word of RM_7(20,4) is planar" (rm20_wW_count.py): random W-point
sets from structured strata of span >= 3, brute-force null spaces of the 35 x W cubic evaluation matrices, and any
full-support null vector on a non-planar set is a violation.
usage: PYTHONPATH=. python scratchpad/rm20_w9_sample_check.py W PER_STRATUM SEED OUT.json [STRATA,comma,separated]"""
import sys, json, time
import numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.pointsets import batched_rank_modp
from qmsd.structured_ad import _null_space

W, per, seed, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
p, m = 7, 4
P = np.asarray(points(m, p), dtype=np.int64)
idx = {tuple(int(v) for v in x): i for i, x in enumerate(P)}
R3 = np.asarray(rm_generator(3, m, p), dtype=np.int64) % p
rng = np.random.default_rng(seed)

def pt(v): return idx[tuple(int(a) % p for a in v)]
def rand_dirs(k):
    while True:
        D = rng.integers(0, p, (k, m))
        if batched_rank_modp(D[None], p)[0] == k: return D
def flat_points(k):
    a = rng.integers(0, p, m); D = rand_dirs(k)
    coeffs = np.array(np.meshgrid(*[np.arange(p)] * k)).reshape(k, -1).T
    return [pt(a + c @ D) for c in coeffs]
def line(a, d): return [pt(a + t * d) for t in range(p)]
def conic_points():
    a = rng.integers(0, p, m); D = rand_dirs(2)
    while True:
        Q = rng.integers(0, p, 6)            # q(s,t) = Q0 s^2 + Q1 st + Q2 t^2 + Q3 s + Q4 t + Q5
        pts = [pt(a + s * D[0] + t * D[1]) for s in range(p) for t in range(p)
               if (Q[0]*s*s + Q[1]*s*t + Q[2]*t*t + Q[3]*s + Q[4]*t + Q[5]) % p == 0]
        if len(pts) >= W - 2 and span(pts) == 2: return pts, a, D     # a genuine plane conic, not a (double) line
def span(S):
    D = (P[S[1:]] - P[S[0]]) % p
    return int(batched_rank_modp(D[None], p)[0])

def stratum(name):
    if name == "uniform":
        return list(rng.choice(p ** m, W, replace=False))
    if name == "3-flat":
        f = flat_points(3); return list(rng.choice(f, W, replace=False))
    if name == "line-heavy":            # 5-7 points of a line + the rest in a random 3-flat through it
        a = rng.integers(0, p, m); D = rand_dirs(3); L = line(a, D[0]); kL = int(rng.integers(5, 8))
        f = [pt(a + c @ D) for c in np.array(np.meshgrid(*[np.arange(p)] * 3)).reshape(3, -1).T]
        rest = [x for x in f if x not in L]
        return list(rng.choice(L, kL, replace=False)) + list(rng.choice(rest, W - kL, replace=False))
    if name == "two-skew-lines":         # points of two skew lines (+ one off point at W = 9 if needed)
        while True:
            a, b = rng.integers(0, p, m), rng.integers(0, p, m); D = rand_dirs(3)
            L1, L2 = line(a, D[0]), line(b, D[1])
            if span(L1 + L2) == 3 and not set(L1) & set(L2): break
        k1 = int(rng.integers(max(1, W - 7), min(7, W - 1) + 1))
        S = list(rng.choice(L1, k1, replace=False)) + list(rng.choice(L2, W - k1, replace=False))
        return S if len(set(S)) == W else stratum(name)
    if name == "conic+off":              # W-1 or W-2 points of a plane conic + points off the plane
        pts, a, D = conic_points(); k = int(rng.integers(W - 2, min(W, len(pts) + 1)))
        if k > len(pts): k = len(pts)
        S = list(rng.choice(pts, k, replace=False))
        while len(S) < W:                # any further points; the caller keeps sets of span >= 3
            x = int(rng.integers(0, p ** m))
            if x not in S: S.append(x)
        return S
    if name == "plane+off":              # W-1 points of a random plane + 1 off it
        f = flat_points(2); S = list(rng.choice(f, W - 1, replace=False))
        while True:
            x = int(rng.integers(0, p ** m))
            if x not in f: return S + [x]
    raise ValueError(name)

T0 = time.time(); res = {}
STRATA = sys.argv[5].split(",") if len(sys.argv) > 5 else ["uniform", "3-flat", "line-heavy", "two-skew-lines", "conic+off", "plane+off"]
for name in STRATA:
    sets = []
    while len(sets) < per:
        S = [int(v) for v in stratum(name)]
        if len(set(S)) == W and span(S) >= 3: sets.append(S)
    sets = np.array(sets)
    rk = batched_rank_modp(R3.T[sets], p)
    dep = sets[rk < W]; full = 0; viol = []
    for S in dep:
        ns = np.asarray(_null_space(R3[:, list(S)], p)) % p
        if ns.shape[0] and ((ns != 0).any(axis=0)).all():
            # no coloop: over F_7 a full-support vector exists
            full += 1; viol.append([tuple(int(v) for v in P[x]) for x in S])
    res[name] = {"sets": int(len(sets)), "dependent": int(len(dep)), "full_support_nonplanar": full, "examples": viol[:5]}
    print(f"[{time.time()-T0:5.0f}s] W={W} {name}: {len(sets)} sets, {len(dep)} dependent, {full} full-support (violations)", flush=True)
json.dump({"W": W, "per_stratum": per, "seed": seed, "strata": res, "violations": sum(v["full_support_nonplanar"] for v in res.values())}, open(out, "w"), indent=1)
