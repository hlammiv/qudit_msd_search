"""Independent random cross-check of the weight-7 lemma for RM_7(20,4) = RM(3,4)^perp.

Draws random 7-subsets of F_7^4 in several strata (uniform on the grid, inside a random plane,
inside a random 3-flat, line-heavy sets, subsets of two lines, subsets of plane conics), computes
the null space of RM(3,4) restricted to each set by brute force (batched F_7 rank of the 35 x 7
matrix, then an explicit null space for the dependent ones) and checks that

  * every full-support null vector lives on a full line ('line7'), and
  * every dependent 7-set contains >= 5 collinear points (the finite Cayley-Bacharach form).

Writes runs/weighted/rm20_w7_sample_check.json.  usage: PYTHONPATH=. python scratchpad/rm20_w7_sample_check.py [N]
"""
import sys, json, time
import numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.pointsets import batched_rank_modp, full_support_null_vector, classify_pointset, affine_lines
from qmsd.structured_ad import _null_space, _rank

p, m = 7, 4
N = int(sys.argv[1]) if len(sys.argv) > 1 else 100_000
T0 = time.time()
def log(s): print(f"[{time.time()-T0:5.0f}s] {s}", flush=True)
P = np.asarray(points(m, p), dtype=np.int64)
POW = p ** np.arange(m)
R3 = np.asarray(rm_generator(3, m, p), dtype=np.int64) % p
rng = np.random.default_rng(20261030)
LINES = affine_lines(p, m)

def choose(n_sets, pool):
    """pool: (n_sets, K) point indices -> random 7-subsets (n_sets, 7)."""
    keys = rng.random(pool.shape)
    order = np.argsort(keys, axis=1)[:, :7]
    return np.take_along_axis(pool, order, axis=1)

def random_flats(n_sets, j):
    """(n_sets, p^j) point indices of random affine j-flats (random base + j independent directions)."""
    out = np.empty((n_sets, p ** j), dtype=np.int64); filled = 0
    Fj = np.asarray(points(j, p), dtype=np.int64)
    while filled < n_sets:
        need = n_sets - filled
        base = rng.integers(0, p, (need, m)); D = rng.integers(0, p, (need, j, m))
        ok = batched_rank_modp(D, p) == j
        pts = (base[ok][:, None, :] + Fj[None, :, :] @ D[ok]) % p         # (ok, p^j, m)
        idx = pts @ POW
        out[filled:filled + len(idx)] = idx; filled += len(idx)
    return out

def stratum_uniform(n):
    return np.array([rng.choice(p ** m, 7, replace=False) for _ in range(n)])

def stratum_plane(n): return choose(n, random_flats(n, 2))
def stratum_3flat(n): return choose(n, random_flats(n, 3))

def stratum_line_heavy(n):
    """s in {4,5,6,7} points of a random line + 7 - s random points off it."""
    out = np.empty((n, 7), dtype=np.int64)
    li = rng.integers(0, len(LINES), n); s = rng.integers(4, 8, n)
    for i in range(n):
        L = LINES[li[i]]; on = rng.choice(L, s[i], replace=False)
        off = []
        while len(off) < 7 - s[i]:
            x = int(rng.integers(0, p ** m))
            if x not in set(L) and x not in off: off.append(x)
        out[i] = np.concatenate([on, np.array(off, dtype=np.int64)])
    return out

def stratum_two_lines(n):
    """7-subsets of L1 u L2 (L2 through a random point, random direction: parallel / meeting / skew)."""
    out = np.empty((n, 7), dtype=np.int64); i = 0
    while i < n:
        L1 = LINES[rng.integers(0, len(LINES))]; L2 = LINES[rng.integers(0, len(LINES))]
        U = np.array(sorted(set(L1.tolist()) | set(L2.tolist())))
        if len(U) < 8: continue
        out[i] = rng.choice(U, 7, replace=False); i += 1
    return out

def stratum_conic(n):
    """7-subsets of the zero set of a random (non-degenerate or not) quadric in a random plane."""
    out = np.empty((n, 7), dtype=np.int64); i = 0; F2 = np.asarray(points(2, p), dtype=np.int64)
    while i < n:
        flat = random_flats(1, 2)[0]                                   # plane points indexed by F2 order
        c = rng.integers(0, p, 6)                                      # 1, s, t, s^2, st, t^2
        s, t = F2[:, 0], F2[:, 1]
        q = (c[0] + c[1] * s + c[2] * t + c[3] * s * s + c[4] * s * t + c[5] * t * t) % p
        Z = flat[q == 0]
        if len(Z) < 7 or len(Z) > 20: continue
        out[i] = rng.choice(Z, 7, replace=False); i += 1
    return out

strata = {"uniform": stratum_uniform, "plane": stratum_plane, "3flat": stratum_3flat,
          "line_heavy": stratum_line_heavy, "two_lines": stratum_two_lines, "conic": stratum_conic}
result = {"N_per_stratum": N, "strata": {}}
for name, fn in strata.items():
    t = time.time()
    dep_labels = {}; full_labels = {}; nullities = {}; n_dep = 0; n_full = 0; bad = []
    for lo in range(0, N, 10_000):                                   # generate and test in chunks (memory)
        T = fn(min(10_000, N - lo))
        assert (np.sort(T, axis=1)[:, 1:] != np.sort(T, axis=1)[:, :-1]).all(), "repeated point"
        rk = batched_rank_modp(R3.T[T], p)
        for Tset in T[rk < 7]:
            n_dep += 1
            ns = np.asarray(_null_space(R3[:, Tset], p)) % p
            nullities[int(ns.shape[0])] = nullities.get(int(ns.shape[0]), 0) + 1
            lab = classify_pointset(P[Tset], p)
            dep_labels[lab] = dep_labels.get(lab, 0) + 1
            u = full_support_null_vector(R3[:, Tset], p)
            if u is not None:
                n_full += 1; full_labels[lab] = full_labels.get(lab, 0) + 1
                assert ((R3[:, Tset] @ u) % p == 0).all() and u.all()
                if lab != "line7": bad.append({"points": P[Tset].tolist(), "label": lab, "word": u.tolist()})
            if not lab.startswith("line"): bad.append({"points": P[Tset].tolist(), "label": lab, "dependent_not_5collinear": True})
    rec = {"sets": N, "dependent": n_dep, "full_support": n_full, "dependent_labels": dep_labels,
           "full_support_labels": full_labels, "nullities": {str(k): v for k, v in nullities.items()},
           "violations": bad, "seconds": round(time.time() - t, 1)}
    result["strata"][name] = rec
    log(f"{name}: {N} sets, dependent {n_dep} {dep_labels}, full-support {n_full} {full_labels}, violations {len(bad)}")
json.dump(result, open("runs/weighted/rm20_w7_sample_check.json", "w"), indent=1)
viol = sum(len(r["violations"]) for r in result["strata"].values())
log(f"RESULT: {viol} violations of the weight-7 lemma in {len(strata) * N} random 7-sets")
assert viol == 0
