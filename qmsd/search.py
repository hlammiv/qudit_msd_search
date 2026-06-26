"""Search for low-overhead qudit triorthogonal codes (NOTES sec 5, sec 10 item 13).

Two engines, matching the paper:
  * manhattan_sweep -- the analytic Manhattan-weight family (exact, integer-only, scales
    to astronomically large blocks); sweep the cutoff w and read off [[n,k,d]] and gamma.
  * random_search    -- randomized puncture-location search over RM_p(r_max,m) (this is how
    the paper found its best small codes, e.g. [[519,106,5]]_5); bounded to small p^m.

search(p) drives both across a range of m and returns the best codes by yield gamma and by
single-round distillation cost C.
"""
from __future__ import annotations

import random

from .reedmuller import r_max, d_rm
from .codes import code_from_manhattan, code_from_puncture, Code
from .distillation import nbar_T, cost

# Largest p^m for which the explicit (distance-certifying) search is run. The
# meet-in-the-middle minimum-distance routine (qmsd.mindist) certifies distances up to 6
# for blocks of this size in seconds, so the explicit search reaches the paper's regime
# (it certified [[519,106,5]]_5 and [[690,39,5]]_3). The analytic Manhattan engine still
# has no size limit. Distances above 6 are left uncertified (the MITM cap).
EXPLICIT_MAX_BLOCK = 750


def manhattan_sweep(p, m, r=None) -> list:
    """Sweep the Manhattan cutoff w; return the valid analytic codes, sorted by gamma.

    Exact and fast (no matrices): for each w it reads n=[m,>w]_p, k=[m,<=w]_p,
    d=Delta_p(m,rtilde,w). Keeps codes with d>=2, k>=1, n>k (so gamma is defined).
    """
    if r is None:
        r = r_max(m, p)
    out: list[Code] = []
    for w in range(0, m * (p - 1)):
        c = code_from_manhattan(p, m, w, r=r)
        if c.d is not None and c.d >= 2 and c.k >= 1 and c.n > c.k:
            out.append(c)
    out.sort(key=lambda c: c.gamma)
    return out


def random_search(p, m, trials, seed=0, target_k=None, max_distance=6) -> list:
    """Randomized puncture-location search over RM_p(r_max,m) (NOTES sec 5).

    Samples ``trials`` puncture-column sets (deterministically, from ``seed``), builds the
    triorthogonal code, and keeps every full-rank, distance-certified candidate (deduped by
    puncture set), sorted by gamma. ``max_distance`` bounds the distance certification so the
    scan stays tractable; codes whose distance exceeds it are skipped, not mis-reported.

    If ``target_k`` is given, every sample has exactly that many punctures; otherwise the
    puncture count is drawn per trial. Bounded to small p**m (it builds F_p matrices).
    """
    rng = random.Random(seed)
    pm = p ** m
    # A loose upper bound on #punctures worth trying (full rank is guaranteed below d_RM,
    # but the check in build_triorthogonal_code is the real gate, so allow a bit more).
    cap = min(pm - 1, max(2, 2 * d_rm(r_max(m, p), m, p)))
    best: dict[frozenset, Code] = {}
    for _ in range(trials):
        k = target_k if target_k is not None else rng.randint(1, cap)
        k = min(k, pm - 1)
        cols = tuple(sorted(rng.sample(range(1, pm + 1), k)))
        key = frozenset(cols)
        if key in best:
            continue
        c = code_from_puncture(p, m, cols, compute_A_d=False, max_distance=max_distance)
        if not c.full_rank or not c.d_certified or c.d is None or c.d < 2 or c.n <= c.k:
            continue
        best[key] = c
    out = list(best.values())
    out.sort(key=lambda c: c.gamma)
    return out


def _cost(c: Code, delta_in: float) -> float:
    """Single-round distillation cost C = n / nbar_T (NOTES eq 39)."""
    return cost(c.n, nbar_T(c.n, c.k, c.p, delta_in))


def search(p, m_values=None, trials_per_m=2000, seed=0, delta_in=1e-3, top=10) -> dict:
    """Top-level: given a prime p, redo the paper's search across a range of m.

    Runs the analytic Manhattan sweep (all m) plus a randomized explicit search (small p**m),
    and returns the best codes ranked by yield gamma and by single-round cost C.
    Returns {best_by_gamma, best_by_cost, scanned}.
    """
    assert isinstance(p, int) and p >= 2
    if m_values is None:
        m_values = [m for m in range(2, 12) if p ** m <= 2000] or [2, 3]

    found: list[Code] = []
    for m in m_values:
        found.extend(manhattan_sweep(p, m))
        # Explicit search only on small blocks: certifying minimum distance scans
        # ~C(p^m, d) column subsets, which is infeasible once p^m grows. The analytic
        # Manhattan engine above covers all sizes (exact, no matrices).
        if p ** m <= EXPLICIT_MAX_BLOCK:
            found.extend(random_search(p, m, trials_per_m, seed=seed))

    # Dedup by (n,k,d) keeping the first occurrence.
    seen: set = set()
    uniq: list[Code] = []
    for c in found:
        key = (c.n, c.k, c.d)
        if c.d is None or key in seen:
            continue
        seen.add(key)
        uniq.append(c)

    by_gamma = sorted(uniq, key=lambda c: c.gamma)[:top]
    by_cost = sorted(uniq, key=lambda c: _cost(c, delta_in))[:top]
    return {
        "p": p,
        "best_by_gamma": by_gamma,
        "best_by_cost": by_cost,
        "scanned": {"m_values": list(m_values), "trials_per_m": trials_per_m,
                    "n_candidates": len(uniq)},
    }
