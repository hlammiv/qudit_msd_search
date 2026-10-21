"""Search for signed weakly self-dual inner codes from affine-flat collections (exact balance -> every Z_{2^k}).

Generators: indicators of affine flats F = {x in F_2^m : A x = b} (codim c = 1..cmax), restricted to the union
of their supports (idle points would be weight-1 logicals).  tau(x) = (-1)^{x_1}.  A flat (or a nonempty
intersection of two flats) is tau-balanced iff e_1 is not in the row space of its equations, so a collection
whose flats and nonempty pairwise intersections all avoid e_1 is signed self-orthogonal over Z.  The signed
punctured RM(r,m) is the collection of coordinate subcubes x_A = 1 (|A| <= r) for coordinates {2..2r+2}; here
any flats are allowed.

Distance (all 2^k): D = min(d_bin(S^perp \\ S), ...) >= min(d_bin(S^perp), wmin(S)) and wmin(S) >= 2^{m-cmax}
(S lies in RM(cmax, m)), so D = d_bin(S^perp) whenever that is < 2^{m-cmax}.  Cost (HHPW, order 5): 1 + 4 n / k.
"""

from __future__ import annotations

import itertools as it
import math

import numpy as np

from qmsd.mindist import min_dependent_columns


def _rank2(rows):
    A = np.array(rows, dtype=np.uint8) % 2
    r = 0
    for c in range(A.shape[1]):
        p = np.nonzero(A[r:, c])[0]
        if not len(p):
            continue
        p = r + p[0]
        A[[r, p]] = A[[p, r]]
        nz = np.nonzero(A[:, c])[0]
        for i in nz:
            if i != r:
                A[i] ^= A[r]
        r += 1
        if r == A.shape[0]:
            break
    return r


class FlatPool:
    def __init__(self, m: int, cmax: int = 2, rng=None, ell=None):
        self.m = m
        self.ell = np.eye(m, dtype=np.uint8)[0] if ell is None else np.asarray(ell, dtype=np.uint8)
        self.pts = np.array([[(v >> j) & 1 for j in range(m)] for v in range(2 ** m)], dtype=np.uint8)
        self.rng = rng if rng is not None else np.random.default_rng()
        self.cmax = cmax

    def random_flat(self, c):
        while True:
            A = self.rng.integers(0, 2, size=(c, self.m)).astype(np.uint8)
            if _rank2(A) < c:
                continue
            if _rank2(np.vstack([A, self.ell[None, :]])) == c:                      # ell in rowspace
                continue
            b = self.rng.integers(0, 2, size=c).astype(np.uint8)
            ind = np.all((self.pts @ A.T) % 2 == b, axis=1)
            return A, b, ind


def compatible(F, G, ell):
    """pairwise tau-balance of F cap G (empty is fine)."""
    A1, b1, ind1 = F; A2, b2, ind2 = G
    if not np.any(ind1 & ind2):
        return True
    A = np.vstack([A1, A2])
    return _rank2(np.vstack([A, np.asarray(ell, dtype=np.uint8)[None, :]])) > _rank2(A)


def subcube(m, A_idx, ell=None):
    pts = np.array([[(v >> j) & 1 for j in range(m)] for v in range(2 ** m)], dtype=np.uint8)
    A = np.zeros((len(A_idx), m), np.uint8)
    for i, a in enumerate(A_idx):
        A[i, a] = 1
    b = np.ones(len(A_idx), np.uint8)
    return A, b, np.all((pts @ A.T) % 2 == b, axis=1)


def signed_rm_flats(r, m):
    """signed punctured RM(r,m) as flats: subcubes x_A = 1, 1 <= |A| <= r; ell = e_1 + ... + e_{2r+1}."""
    ell = np.zeros(m, np.uint8); ell[:2 * r + 1] = 1
    return [subcube(m, A) for deg in range(1, r + 1) for A in it.combinations(range(m), deg)], ell


def evaluate(flats, m, D_max=6, ell=None):
    ind = np.array([f[2] for f in flats], dtype=np.uint8)
    support = np.nonzero(ind.any(0))[0]
    S = ind[:, support]
    n = len(support)
    r = _rank2(S)
    k = n - 2 * r
    pts = np.array([[(v >> j) & 1 for j in range(m)] for v in support], dtype=np.int64)
    ell = np.eye(m, dtype=np.int64)[0] if ell is None else np.asarray(ell, dtype=np.int64)
    tau = 1 - 2 * ((pts @ ell) % 2)
    B = (S.astype(np.int64) * tau) @ S.T.astype(np.int64)
    assert not np.any(B), "unbalanced"
    if k <= 0:
        return {"n": n, "k": k, "D": None}
    D_max = min(D_max, 6)
    try:
        d = min_dependent_columns(S.astype(np.int64), 2, d_max=D_max)
    except ValueError:
        d = None
    D = d if d is not None and d <= D_max else D_max + 1
    odd = None
    normal = _rank2(list(S) + [np.ones(n, np.uint8)]) > r          # 1 not in S  <=>  an odd-weight logical exists
    return {"n": int(n), "k": int(k), "r": int(r), "D": int(D), "D_exact": D <= D_max, "normal": bool(normal),
            "cost5": 1 + 4 * n / k, "support": support.tolist(),
            "rows": ["".join(map(str, row)) for row in S], "ell": ell.tolist()}


def greedy(m, target_rank, cmax=2, rng=None, max_tries=4000, codim_p=None, ell=None):
    pool = FlatPool(m, cmax, rng, ell)
    rng = pool.rng
    flats, rows = [], []
    tries = 0
    while len(flats) < target_rank and tries < max_tries:
        tries += 1
        c = int(rng.choice(np.arange(1, cmax + 1), p=codim_p))
        F = pool.random_flat(c)
        if not all(compatible(F, G, pool.ell) for G in flats):
            continue
        if _rank2(rows + [F[2]]) <= len(rows):
            continue
        flats.append(F); rows.append(F[2])
    return flats


def anneal(m, start, ell, steps=3000, k_range=(16, 24), cmax=2, seed=0, T0=0.5, D_min=5, log=None, s_out=2):
    """Simulated annealing over flat collections: moves = drop a flat / add a compatible random flat / swap.
    Score = cost5 (1 + 4n/k) for D >= D_min and k in k_range, else penalised."""
    rng = np.random.default_rng(seed)
    pool = FlatPool(m, cmax, rng, ell)

    def score(fl):
        if not fl:
            return 1e9, None
        rows = [f[2] for f in fl]
        if _rank2(rows) < len(rows):
            return 1e9, None
        e = evaluate(fl, m, min(D_min, 6), ell)
        if e["k"] <= 0:
            return 1e9, e
        pen = 0.0
        if e["D"] < D_min:
            pen += 50 * (D_min - e["D"])
        if not e["normal"]:
            pen += 100
        lo, hi = k_range
        if e["k"] < lo:
            pen += 2 * (lo - e["k"])
        if e["k"] > hi:
            pen += 2 * (e["k"] - hi)
        return 1 + 2 * s_out * e["n"] / e["k"] + pen, e

    cur = list(start); cs, ce = score(cur)
    best = (cs, ce, list(cur))
    for step in range(steps):
        T = T0 * (1 - step / steps) + 1e-3
        new = list(cur)
        mv = rng.random()
        if mv < 0.4 and len(new) > 1:
            new.pop(int(rng.integers(len(new))))
        else:
            if mv > 0.7 and len(new) > 1:
                new.pop(int(rng.integers(len(new))))
            for _ in range(50):
                F = pool.random_flat(int(rng.integers(1, cmax + 1)))
                if all(compatible(F, G, ell) for G in new):
                    new.append(F); break
        ns, ne = score(new)
        if ns < cs or rng.random() < math.exp(-(ns - cs) / T):
            cur, cs, ce = new, ns, ne
            if cs < best[0]:
                best = (cs, ce, list(cur))
                if log:
                    log(step, ce)
    return best
