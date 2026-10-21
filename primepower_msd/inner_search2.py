"""Inner-code search v2: annealing over signed flat collections with a graded distance objective.

The v1 annealer (inner_search.anneal) only sees "binary d(S^perp) <= 6 or not", so any flat removal from
signed RM(2,m) is a cliff.  Here the objective counts N_bad = #{pairs of distinct column subsets of size <= 3
with equal syndrome} (empty set included): N_bad = 0  <=>  no nonzero dual word of weight <= 6  <=>
d_bin(S^perp) >= 7, and it decreases smoothly as low-weight dual words are removed.

Generators: affine flats of codim 1..cmax in F_2^m, pairwise tau-compatible for tau = (-1)^{ell.x}
(inner_search.compatible), restricted to the union of their supports.  Exact integer balance is asserted.

    PYTHONPATH=. python -m primepower_msd.inner_search2 --m 7 --kmin 35 --kmax 50 --seed 0 --out runs/meas_distill/inner7b/a_s0.json
"""
from __future__ import annotations

import argparse
import json
import math
import time

import numpy as np
from numba import njit

from .inner_search import FlatPool, compatible, signed_rm_flats, _rank2


@njit(cache=False)
def _basis_syndromes(S):
    """Row-reduce 0/1 matrix S (R x n) over F_2; return (rank, column syndromes as uint64 (rank <= 64),
    ones_in_span)."""
    A = S.copy()
    R, n = A.shape
    r = 0
    piv = np.empty(R, np.int64)
    for c in range(n):
        p = -1
        for i in range(r, R):
            if A[i, c]:
                p = i; break
        if p < 0:
            continue
        if p != r:
            for j in range(n):
                t = A[r, j]; A[r, j] = A[p, j]; A[p, j] = t
        for i in range(R):
            if i != r and A[i, c]:
                for j in range(n):
                    A[i, j] ^= A[r, j]
        piv[r] = c
        r += 1
        if r == R:
            break
    syn = np.zeros(n, np.uint64)
    rr = min(r, 64)
    for i in range(rr):
        for j in range(n):
            if A[i, j]:
                syn[j] |= np.uint64(1) << np.uint64(i)
    # all-ones in row space: reduce ones by pivots
    v = np.ones(n, np.uint8)
    for i in range(r):
        if v[piv[i]]:
            for j in range(n):
                v[j] ^= A[i, j]
    ones_in = True
    for j in range(n):
        if v[j]:
            ones_in = False; break
    return r, syn, ones_in


@njit(cache=False)
def _nbad(syn):
    n = syn.shape[0]
    N = 1 + n + n * (n - 1) // 2 + n * (n - 1) * (n - 2) // 6
    h = np.empty(N, np.uint64)
    t = 0
    h[t] = 0; t += 1
    for i in range(n):
        h[t] = syn[i]; t += 1
    for i in range(n):
        for j in range(i + 1, n):
            h[t] = syn[i] ^ syn[j]; t += 1
    for i in range(n):
        for j in range(i + 1, n):
            a = syn[i] ^ syn[j]
            for k in range(j + 1, n):
                h[t] = a ^ syn[k]; t += 1
    h.sort()
    tot = 0
    run = 1
    for i in range(1, N):
        if h[i] == h[i - 1]:
            run += 1
        else:
            tot += run * (run - 1) // 2
            run = 1
    tot += run * (run - 1) // 2
    return tot


def evaluate(flats, m, ell):
    ind = np.array([f[2] for f in flats], dtype=np.uint8)
    support = np.nonzero(ind.any(0))[0]
    S = np.ascontiguousarray(ind[:, support])
    n = len(support)
    pts = np.array([[(v >> j) & 1 for j in range(m)] for v in support], dtype=np.int64)
    tau = 1 - 2 * ((pts @ np.asarray(ell, np.int64)) % 2)
    B = (S.astype(np.int64) * tau) @ S.T.astype(np.int64)
    assert not np.any(B), "unbalanced"
    r, syn, ones_in = _basis_syndromes(S)
    k = n - 2 * r
    e = {"n": int(n), "k": int(k), "r": int(r), "normal": not ones_in, "support": support, "S": S}
    e["nbad"] = int(_nbad(syn)) if r <= 64 else 10 ** 9
    return e


def to_json(e, ell):
    return {"n": e["n"], "k": e["k"], "r": e["r"], "nbad": e["nbad"], "D_bin_ge7": e["nbad"] == 0,
            "normal": e["normal"], "n_over_k": e["n"] / e["k"], "T_per_out7": 1 + 6 * e["n"] / e["k"],
            "gamma7": math.log(1 + 6 * e["n"] / e["k"]) / math.log(7),
            "support": e["support"].tolist(), "rows": ["".join(map(str, row)) for row in e["S"]],
            "ell": [int(x) for x in ell]}


def score(e, k_range, alpha):
    if e["k"] <= 0:
        return 1e9
    lo, hi = k_range
    pen = alpha * math.sqrt(e["nbad"]) + (0 if e["normal"] else 20.0)
    pen += 1.0 * max(0, lo - e["k"]) + 1.0 * max(0, e["k"] - hi)
    return 1 + 6 * e["n"] / e["k"] + pen


def anneal(m, start, ell, steps, k_range, cmax, rng, T0=1.0, alpha=0.5, record=None, cprobs=None, forbid=None):
    """forbid: boolean mask over F_2^m of points no generator may touch (caps n at 2^m - |forbid|)."""
    pool = FlatPool(m, cmax, rng, ell)
    if forbid is not None:
        start = [F for F in start if not np.any(F[2] & forbid)]
        if not start:
            for _ in range(100000):
                F = pool.random_flat(int(rng.integers(1, cmax + 1)))
                if F[2].any() and not np.any(F[2] & forbid):
                    start = [F]; break
    cs_ = np.arange(1, cmax + 1)
    cprobs = np.ones(cmax) / cmax if cprobs is None else np.asarray(cprobs) / np.sum(cprobs)
    cur = list(start)
    ce = evaluate(cur, m, ell); cs = score(ce, k_range, alpha)
    best = (cs, ce)
    for step in range(steps):
        T = T0 * (1 - step / steps) + 1e-3
        new = list(cur)
        mv = rng.random()
        if mv < 0.3 and len(new) > 2:
            new.pop(int(rng.integers(len(new))))
        else:
            if mv > 0.7 and len(new) > 2:
                new.pop(int(rng.integers(len(new))))
            for _ in range(100):
                F = pool.random_flat(int(rng.choice(cs_, p=cprobs)))
                if not F[2].any() or (forbid is not None and np.any(F[2] & forbid)):
                    continue
                if all(compatible(F, G, ell) for G in new):
                    new.append(F); break
        if _rank2([f[2] for f in new]) < len(new):
            continue
        ne = evaluate(new, m, ell)
        ns = score(ne, k_range, alpha)
        if ns < cs or rng.random() < math.exp(-(ns - cs) / T):
            cur, cs, ce = new, ns, ne
            if record:
                record(step, ce)
            if cs < best[0]:
                best = (cs, ce)
    return best


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, default=7)
    ap.add_argument("--kmin", type=int, default=35)
    ap.add_argument("--kmax", type=int, default=50)
    ap.add_argument("--cmax", type=int, default=3)
    ap.add_argument("--cprobs", type=float, nargs="+", default=None)
    ap.add_argument("--steps", type=int, default=5000)
    ap.add_argument("--restarts", type=int, default=4)
    ap.add_argument("--start", choices=["rm", "thin", "empty"], default="rm")
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--T0", type=float, default=1.0)
    ap.add_argument("--forbid", type=int, default=0, help="codim of a random affine flat of forbidden points (0 = none)")
    ap.add_argument("--aligned", action="store_true", help="forbidden flat = {x_j = 0 for the last `forbid` coords}")
    ap.add_argument("--forbid_extra", type=int, default=0, help="extra random forbidden points")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    rng = np.random.default_rng(a.seed)
    found = {}
    t0 = time.time()

    def record(step, e):
        if e["nbad"] == 0 and e["normal"] and a.kmin <= e["k"] <= a.kmax:
            k = e["k"]
            if k not in found or e["n"] < found[k]["n"]:
                found[k] = to_json(e, ell)
                print(f"k {k} n {e['n']} r {e['r']} T/out {found[k]['T_per_out7']:.3f} "
                      f"gamma {found[k]['gamma7']:.4f}  step {step} t {time.time() - t0:.0f}s", flush=True)
                json.dump({"m": a.m, "seed": a.seed, "args": vars(a), "found": found}, open(a.out, "w"))

    for rs in range(a.restarts):
        start, ell = signed_rm_flats(2, a.m) if a.m >= 5 else (None, None)
        ell = np.asarray(ell, np.uint8)
        if a.start == "thin":
            idx = rng.permutation(len(start))[: int(rng.integers(len(start) // 2, len(start)))]
            start = [start[i] for i in idx]
        elif a.start == "empty":
            start = start[:1]
        forbid = None
        if a.forbid or a.forbid_extra:
            pts = np.array([[(v >> j) & 1 for j in range(a.m)] for v in range(2 ** a.m)], dtype=np.uint8)
            forbid = np.zeros(2 ** a.m, bool)
            if a.forbid and a.aligned:
                forbid |= np.all(pts[:, a.m - a.forbid:] == 0, axis=1)
            elif a.forbid:
                while True:
                    A = rng.integers(0, 2, size=(a.forbid, a.m)).astype(np.uint8)
                    if _rank2(A) == a.forbid:
                        break
                forbid |= np.all((pts @ A.T) % 2 == rng.integers(0, 2, size=a.forbid), axis=1)
            if a.forbid_extra:
                forbid[rng.choice(np.nonzero(~forbid)[0], a.forbid_extra, replace=False)] = True
        best = anneal(a.m, start, ell, a.steps, (a.kmin, a.kmax), a.cmax, rng, T0=a.T0, alpha=a.alpha,
                      record=record, cprobs=a.cprobs, forbid=forbid)
        print("restart", rs, "best score", round(best[0], 3), {kk: best[1][kk] for kk in ("n", "k", "r", "nbad", "normal")},
              flush=True)
        json.dump({"m": a.m, "seed": a.seed, "args": vars(a), "found": found}, open(a.out, "w"))


if __name__ == "__main__":
    main()
