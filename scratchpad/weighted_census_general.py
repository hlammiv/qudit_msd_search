"""Weighted-support census for a prime p and m: for each r with a non-negative charge budget and
each diagonal quadric / low-degree charge within budget, each cube class: |X|, class points N_c,
class rank rho (k ceiling), max class points on a plane, and the random-count model's best gamma:
for w = 4..8, k_w = largest k with E[accidental block words of weight <= w]
= sum_{v<=w} C(N_c-k, v) (p-1)^(v-1) / p^(rho-k) < 0.1, and gamma_w = ln((|X|-k_w)/k_w)/ln(w+1)
when (w+1) > (|X|-k_w)/k_w. Reports the best (w, k_w, gamma_w) per row.   usage: p m [maxdeg]"""
import sys, itertools, numpy as np
from math import comb, log
from qmsd.reedmuller import points
from qmsd.weighted import weighted_support, quadric, class_points, charge_degree_budget, product_of_linear
from qmsd.structured_ad import _rank
p, m = int(sys.argv[1]), int(sys.argv[2]); maxdeg = int(sys.argv[3]) if len(sys.argv) > 3 else 2
P = np.asarray(points(m, p), dtype=np.int64)
units = sorted({pow(u, 3, p) and u for u in range(1, p)})  # representatives: 1 plus non-cubes
cubes = {pow(x, 3, p) for x in range(1, p)}
reps = [1] + [u for u in range(2, p) if u not in cubes][: (p - 1) // len(cubes) - 1]
def model_best(Nc, rho, nX, rk):
    """Accidental words live either in the class block (N_c - k columns, rho - k checks) or anywhere
    (|X| - k columns, rank - k checks); both populations are counted."""
    best = None
    for w in range(4, 9):
        kw = None
        for k in range(10, rho):
            E = sum((comb(Nc - k, v) / p ** (rho - k) + comb(nX - k, v) / p ** (rk - k)) * (p - 1) ** (v - 1) for v in range(2, w + 1))
            if E < 0.1: kw = k
        if kw is None: continue
        nk = (nX - kw) / kw
        if w + 1 > nk:
            g = log(nk) / log(w + 1)
            if best is None or g < best[2]: best = (w + 1, kw, g)
    return best
def max_on_plane(idx):
    """max number of the given points (indices into the grid) on an affine 2-flat (m>=3) or line (m=2)."""
    pts = P[idx]; best = 0
    dirs = [v for v in itertools.product(range(p), repeat=m) if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1]
    if m == 2:
        for a in dirs:
            vals = (pts @ np.array(a)) % p; best = max(best, int(np.bincount(vals, minlength=p).max()))
        return best
    seen = set()
    for a in dirs:
        for b in dirs:
            if a >= b: continue
            key = (a, b)
            M = np.array([a, b]); if_rank = _rank(M.copy(), p)
            if if_rank < 2: continue
            v = (pts @ M.T) % p; k = v[:, 0] * p + v[:, 1]; best = max(best, int(np.bincount(k, minlength=p * p).max()))
    return best
rmax = (m * (p - 1) - 1) // 3
rows = []
print(f"p={p} m={m}: grid {p**m}; r_max(unweighted)={rmax}", flush=True)
for r in range(rmax, max(rmax - 3, 0), -1):
    B = charge_degree_budget(p, m, r)
    if B < 0: continue
    forms = []
    if B >= 1:
        for lin in itertools.product(range(2), repeat=m):
            if any(lin):
                for c in range(p): forms.append((f"lin{lin}+{c}", (P @ np.array(lin) + c) % p))
    if B >= 2:
        nonsq = next(u for u in range(2, p) if pow(u, (p - 1) // 2, p) != 1) if p > 2 else 1
        diags = sorted({tuple(sorted(d, reverse=True)) for d in itertools.product((0, 1, nonsq), repeat=m) if any(d)})
        for d in diags:
            for c in range(p): forms.append((f"q{d}+{c}", quadric(p, m, d, const=c)))
            if 0 in d: forms.append((f"q{d}+x", quadric(p, m, d, lin=tuple(1 if x == 0 else 0 for x in d))))
    if B >= 3 and maxdeg >= 3:
        for c in range(p): forms.append((f"cubic sum x^3+{c}", (np.sum(P ** 3, axis=1) + c) % p))
    print(f"-- r={r}: budget deg g <= {B}, {len(forms)} charges", flush=True)
    for name, g in forms:
        if not g.any(): continue
        GX, X, mu = weighted_support(p, m, g, r); rk = _rank(GX, p)
        for u in reps:
            cls = class_points(mu, p, u)
            if len(cls) < 20: continue
            rho = _rank(GX[:, cls], p); mb = model_best(len(cls), rho, len(X), rk)
            mp = max_on_plane([int(X[i]) for i in cls]) if (m <= 4 and len(cls) <= 3000) else -1
            tag = f"d={mb[0]} k={mb[1]} gamma={mb[2]:.3f}" if mb else "no gamma<1 in model"
            print(f"r={r} {name:28s} |X|={len(X):5d} rank={rk:4d} unit={u} N_c={len(cls):5d} rho={rho:4d} vanish={rk-rho:3d} max_codim2_flat={mp:3d} | {tag}", flush=True)
