"""Step 2: fault-tolerant (across-block) ring-CCZ codes — the cheap resource of the two-tier plan.

Setting.  Three copies of a CSS Z_4 code lifted from binary data (S = 0/1 X-stabilizer rows, L = 0/1
logical-X representatives, M_Z = (S + L)^perp over Z_4), with the physical gate
    G = prod_i CCZ_ring(x_i, y_i, z_i),   CCZ_ring = diag(i^{xyz})
between corresponding qudits of the three blocks (every qudit in exactly one 3-body term with legs in
distinct blocks: manifestly fault tolerant, output error ~ A_d p^d).

Hosting (derived; multilinear over Z_4 so generators suffice; confirmed by rgt.hosts):
    Delta_{(a,0,0)} Phi = 4 sum_i a_i y_i z_i == 0 (mod 16)  for a in S, y,z in <S + L>
    <=>  |a ^ u ^ v| == 0 (mod 4)   for all a in S and u, v in S u L      ..... (C4)
Induced logical on the 3k logical ququarts:  F(p,q,r) = sum_{j,j',j''} p_j q_j' r_j'' |l_j ^ l_j' ^ l_j''|
(mod 4): a sum of logical ring-CCZ / CS-type terms whose coefficients are the triple overlaps of the
logical representatives; the gate is a GENUINE order-4 ring-CCZ iff some |l_j ^ l_j' ^ l_j''| is odd.
(Compare qubit CCZ codes: Haah-Hastings need the same overlaps mod 2; the ring needs mod 4.)

Search.  For a binary CSS code (S, Z-checks) with k logicals, the representatives l_j are free to move
by stabilizers (l_j -> l_j + s) and logicals may be DROPPED (their Z partner joins M_Z; k decreases).
`optimize_reps` does a randomized local search over representative choices and subsets to satisfy (C4)
with as many logicals as possible and an odd logical triple overlap; the result is lifted with
lifts.lift_css_from_x, hosting re-verified with rgt, distance certified with z4_distance, gamma recorded.

Families: CSS Hamming [[2^r-1, 2^r-1-2r, 3]], punctured RM(1,m) [[2^m-1,1,3]], Haah-Hastings-type RM(r,m)
(3r < m) puncture codes, random binary CSS codes, and anything in lifts.lift_families.
"""

from __future__ import annotations

import itertools as it
import math
import random
import time

import numpy as np

from . import rgt
from .lifts import binary_rank, lift_css_from_x, rm_binary, rm_puncture_candidate, _binary_dual
from .symplectic import Z4Code, blocks
from .z4_distance import code_distance


# --------------------------------------------------------------------------- the C4 condition
def c4_violations(S: np.ndarray, L: np.ndarray) -> int:
    """Number of (a, u, v) with a in S, u <= v in S u L and |a ^ u ^ v| != 0 mod 4."""
    S = np.atleast_2d(S) % 2
    U = np.vstack([S, np.atleast_2d(L) % 2]) if len(L) else S
    bad = 0
    for a in S:
        M = (U * a[None, :])                     # rows u restricted to supp(a)
        T = M @ M.T                              # |a ^ u ^ v|
        bad += int(np.count_nonzero(np.triu(T % 4)))
    return bad


def logical_triple_odd(L: np.ndarray) -> bool:
    L = np.atleast_2d(L) % 2
    k = len(L)
    for i in range(k):
        for j in range(i, k):
            for l in range(j, k):
                if int(np.sum(L[i] & L[j] & L[l])) % 2:
                    return True
    return False


# --------------------------------------------------------------------------- representative search
def optimize_reps(S: np.ndarray, L0: np.ndarray, rng: random.Random, iters: int = 4000,
                  min_k: int = 1) -> tuple[np.ndarray, int]:
    """Randomized local search: move representatives by stabilizers and drop logicals until (C4) holds
    with an odd logical triple overlap.  Returns (L, violations); violations == 0 means success."""
    S = np.atleast_2d(S) % 2
    L = np.atleast_2d(L0) % 2
    best = (c4_violations(S, L), L.copy())
    for _ in range(iters):
        if best[0] == 0 and logical_triple_odd(best[1]):
            break
        L = best[1].copy()
        move = rng.random()
        if move < 0.7 and len(L):
            j = rng.randrange(len(L))
            s = S[rng.randrange(len(S))]
            L[j] = (L[j] + s) % 2
        elif len(L) > min_k:
            j = rng.randrange(len(L))
            L = np.delete(L, j, axis=0)
        else:
            continue
        v = c4_violations(S, L)
        score = v if not (v == 0 and not logical_triple_odd(L)) else 0.5
        if score <= best[0] + (0 if rng.random() < 0.9 else 1):
            best = (score, L)
    return best[1], best[0]


# --------------------------------------------------------------------------- binary sources
def hamming_css(r: int) -> dict:
    """CSS Hamming [[2^r - 1, 2^r - 1 - 2r, 3]]: X- and Z-checks both the Hamming parity-check matrix."""
    n = 2 ** r - 1
    H = np.array([[(c >> i) & 1 for c in range(1, n + 1)] for i in range(r)], dtype=np.int64)
    Zb = _binary_dual(H)                      # all of H^perp ...
    # logical X reps: a basis of H^perp modulo <H>
    G = np.vstack([H, Zb])
    reps = []
    span_rows = [h for h in H]
    for v in Zb:
        if binary_rank(np.array(span_rows + [v])) > binary_rank(np.array(span_rows)):
            reps.append(v); span_rows.append(v)
    return {"n": n, "k": len(reps), "S": H, "L": np.array(reps, dtype=np.int64), "label": f"Hamming CSS [[{n},{len(reps)},3]]"}


def random_binary_css(n: int, r: int, rng: random.Random, d_min: int = 3, tries: int = 200) -> dict | None:
    """Random binary CSS code with X-checks = Z-checks = H (r x n), requiring d >= d_min via qmsd MITM."""
    from qmsd.mindist import min_dependent_columns
    for _ in range(tries):
        H = np.array([[rng.randrange(2) for _ in range(n)] for _ in range(r)], dtype=np.int64)
        if binary_rank(H) < r or np.any(~np.any(H, axis=0)):
            continue
        if np.any((H @ H.T) % 2):
            continue
        try:
            d = min_dependent_columns(H, 2, d_max=d_min)
        except ValueError:
            d = d_min + 1
        if d is None or d < d_min:
            continue
        Zb = _binary_dual(H)
        reps = []
        span_rows = [h for h in H]
        for v in Zb:
            if binary_rank(np.array(span_rows + [v])) > binary_rank(np.array(span_rows)):
                reps.append(v); span_rows.append(v)
        return {"n": n, "k": len(reps), "S": H, "L": np.array(reps, dtype=np.int64), "label": f"random CSS [[{n},{len(reps)},>={d_min}]]"}
    return None


def rm_puncture_source(m: int, r: int, punctures) -> dict | None:
    c = rm_puncture_candidate(m, punctures, r)
    if c is None:
        return None
    return {"n": c["n"], "k": c["k"], "S": c["MX"], "L": c["Xlog"], "label": f"RM({r},{m}) punct k={c['k']}"}


# --------------------------------------------------------------------------- evaluation
def evaluate_ccz(src: dict, rng: random.Random, iters: int = 4000, d_max: int = 4) -> dict:
    t0 = time.time()
    S, L0 = np.atleast_2d(src["S"]), np.atleast_2d(src["L"])
    L, viol = optimize_reps(S, L0, rng, iters=iters)
    out = {"label": src["label"], "n": src["n"], "k_binary": src["k"], "k_kept": int(len(L)), "c4_violations": int(viol),
           "odd_triple": bool(logical_triple_odd(L)) if len(L) else False}
    if viol != 0 or not len(L):
        out["seconds"] = round(time.time() - t0, 1)
        return out
    code = lift_css_from_x(S, L, src["n"], label=src["label"] + f" k={len(L)}")
    out.update({"k": code.k, "k_genuine": code.k_genuine, "soft": code.soft, "free": code.free_quotient_ok})
    tri = blocks(code, 3)
    layer = rgt.TransversalLayer.across_block_ccz(code.n)
    h = rgt.hosts(tri, layer, "support")
    out["hosted"] = bool(h.hosted)
    out["protocol"] = layer.degree_report([0] * code.n + [1] * code.n + [2] * code.n)["protocol"]
    try:
        dmin, cert, det = code_distance(code, d_max=d_max)
    except MemoryError:
        dmin, cert, det = None, False, {}
    out.update({"distance": dmin, "distance_certified": cert, "distance_detail": det})
    out["gamma"] = (math.log(code.n / code.k) / math.log(dmin)) if (dmin and dmin > 1 and code.k > 0) else None
    if h.hosted and tri.k_genuine <= 6:
        rep = rgt.induced_logical_report(tri, layer)
        out["logical"] = {k: rep[k] for k in ("level", "odd_third_difference", "escapes_A2", "entangling")}
        out["logical"]["ccz_terms"] = len(rep["ccz_terms"])
    out["L"] = L.tolist()
    out["seconds"] = round(time.time() - t0, 1)
    return out


def run_families(seed: int = 0, d_max: int = 4, verbose: bool = True) -> list[dict]:
    rng = random.Random(seed)
    srcs = []
    for r in (3, 4, 5):
        srcs.append(hamming_css(r))
    for m in (4, 5):
        from .lifts import punctured_rm_triorthogonal
        b = punctured_rm_triorthogonal(m)
        srcs.append({"n": b["n"], "k": 1, "S": b["MX"], "L": np.ones((1, b["n"]), dtype=np.int64), "label": b["label"]})
    for (m, r, k) in ((5, 1, 2), (6, 2, 4), (6, 2, 8), (7, 2, 10)):
        for _ in range(3):
            P = sorted(rng.sample(range(2 ** m), k))
            s = rm_puncture_source(m, r, P)
            if s:
                srcs.append(s)
    for (n, r) in ((12, 3), (16, 4), (20, 5), (24, 6)):
        s = random_binary_css(n, r, rng)
        if s:
            srcs.append(s)
    results = []
    for src in srcs:
        res = evaluate_ccz(src, rng, d_max=d_max)
        results.append(res)
        if verbose:
            print({k: v for k, v in res.items() if k not in ("L",)}, flush=True)
    return results


if __name__ == "__main__":
    import json, sys
    res = run_families()
    with open("runs/diverse/catalog_ccz.json", "w") as f:
        json.dump(res, f, indent=1, default=str)
    good = [r for r in res if r.get("hosted") and r.get("odd_triple") and (r.get("distance") or 0) >= 2]
    print("FT ring-CCZ codes (hosted, odd triple, d>=2):")
    for r in sorted(good, key=lambda r: r["gamma"] or 99):
        print(f"  {r['label']:35s} [[{r['n']},{r['k']},{r['distance']}]] gamma={r['gamma']}")
