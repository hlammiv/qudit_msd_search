"""Search for low-overhead qudit triorthogonal codes (NOTES sec 5, sec 10 item 13).

Two engines, matching the paper:
  * manhattan_sweep -- the analytic Manhattan-weight family (exact, integer-only, scales
    to astronomically large blocks); sweep the cutoff w and read off [[n,k,d]] and gamma.
  * random_search    -- randomized puncture-location search over RM_p(r_max,m). Supports
    process-level parallelism (``n_jobs``, trials are independent) and a ``sampler`` mode:
      - "uniform"      : i.i.d. random puncture sets (stalls at d=2 on hard cases).
      - "capset"       : draw cap sets (no 3 collinear points) -- the structure that the
                         paper's high-distance codes have (see SAMPLING_INVESTIGATION.md).
      - "capset_climb" : a cap seed followed by a cap-preserving distance hill-climb; the
                         most efficient at reaching the rare d>=3 puncture sets.

search(p) drives both across a range of m and returns the best codes by yield gamma and by
single-round distillation cost C.
"""
from __future__ import annotations

import os
import random
from dataclasses import replace

from .reedmuller import r_max, d_rm, rm_generator
from .codes import code_from_manhattan, code_from_puncture, Code
from .distillation import nbar_T, cost
from .sampling import all_points, random_cap, cap_extends, points_to_columns, random_plane_spread

EXPLICIT_MAX_BLOCK = 750  # largest p^m for which the distance-certifying explicit search runs
SAMPLERS = ("uniform", "capset", "capset_climb", "arc_climb", "plane_spread")
# arc_climb ranks candidates by the exact A_d surrogate; it routes A_d through the exact
# MacWilliams engine, feasible only when p**dim(G0) <= this budget (the small-dual regime).
_ARC_EXACT_BUDGET = 5_000_000


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


def _search_chunk(p, m, r, cap, pm, n_trials, seed, target_k, max_distance,
                  sampler, climb_steps, swap_tries) -> dict:
    """Evaluate ``n_trials`` candidate puncture sets; return {frozenset(cols): Code} of valid
    (full-rank, distance-certified, d>=2) codes. Self-contained for parallel workers (builds
    its own RM generator and RNG). ``sampler`` selects how each candidate is drawn."""
    rng = random.Random(seed)
    G = rm_generator(r, m, p)
    best: dict = {}

    if sampler == "uniform":
        for _ in range(n_trials):
            k = target_k if target_k is not None else rng.randint(1, cap)
            k = min(k, pm - 1)
            cols = tuple(sorted(rng.sample(range(1, pm + 1), k)))
            key = frozenset(cols)
            if key in best:
                continue
            c = code_from_puncture(p, m, cols, r=r, compute_A_d=False,
                                   max_distance=max_distance, G=G)
            if not c.full_rank or not c.d_certified or c.d is None or c.d < 2 or c.n <= c.k:
                continue
            best[key] = c
        return best

    # --- structure-aware cap-set samplers ---
    allpts = all_points(m, p)
    # arc_climb ranks candidates by the A_d surrogate and so needs A_d; it routes A_d through
    # the exact MacWilliams engine (fast + distance-uncapped when the dual is small). The
    # other structure-aware samplers leave A_d off (cheaper) and rank on distance alone.
    want_ad = sampler == "arc_climb"
    ad_budget = _ARC_EXACT_BUDGET if want_ad else 0

    def _eval(cols):
        """Build/keep the code; return (d_int, full_rank, A_d_or_None, Code)."""
        key = frozenset(cols)
        cached = best.get(key)
        if cached is not None:
            return cached.d, True, cached.A_d, cached
        c = code_from_puncture(p, m, cols, r=r, compute_A_d=want_ad,
                               max_distance=max_distance, G=G, exact_budget=ad_budget)
        ad = c.A_d
        # arc_climb surrogate fallback: when the exact MacWilliams A_d declines (large dual,
        # e.g. m>=5) count the minimum-weight dual codewords DIRECTLY by meet-in-the-middle.
        # That count is dim(G0)-independent, so the (distance, -A_d) gradient survives exactly
        # where MacWilliams dies -- this is what lets arc_climb push d=5 -> d=6 at m=5.
        if want_ad and ad is None and c.d is not None and c.d_certified and 2 <= c.d <= 6:
            from .weightcount import count_weight_d
            from .triorthogonal import build_triorthogonal_code
            try:
                G0 = build_triorthogonal_code(p, m, r, cols, G=G)["X_stab"]
                ad = count_weight_d(G0, c.d, p)
            except (NotImplementedError, MemoryError):
                ad = None
        if c.full_rank and c.d_certified and c.d is not None and c.d >= 2 and c.n > c.k:
            if ad is not None and c.A_d is None:
                c = replace(c, A_d=ad)   # attach the counted A_d to the stored Code
            best[key] = c
        return (c.d if c.d is not None else 0), bool(c.full_rank), ad, c

    # Climb fitness: maximise the certified distance, then MINIMISE the multiplicity of the
    # minimum-weight dual codewords A_d. The integer distance is a flat plateau with isolated
    # d+1 spikes, so on "arc_climb" the A_d term supplies the gradient that carries a cap down
    # its A_d basin until the minimum-weight codewords vanish and the distance ticks up. When
    # A_d is unavailable (None -- dual too large for the exact engine) the tuple degrades to
    # distance-only, i.e. "capset_climb" behaviour. "capset_climb" never computes A_d, so it
    # always ranks on distance alone -- its original behaviour, bit-for-bit.
    HUGE = pm * pm + 1

    def _fit(d_int, full_rank, ad):
        if not full_rank:
            return (-1, 0)
        return (d_int, -(ad if ad is not None else HUGE))

    for _ in range(n_trials):
        k = target_k if target_k is not None else rng.randint(1, cap)
        k = min(k, pm - 1)
        if sampler == "plane_spread":     # cap + no-4-coplanar: reaches the higher-distance codes
            seed_pts = random_plane_spread(m, p, k, rng, allpts)
        else:
            seed_pts = random_cap(m, p, k, rng, allpts)
        if seed_pts is None:  # greedy pass stalled (k too large for the structure); retry next trial
            continue
        cur_pts = seed_pts
        d0, fr0, ad0, _ = _eval(points_to_columns(cur_pts, p))
        if sampler in ("capset", "plane_spread"):  # seed-only samplers: no climb
            continue
        # cap-preserving, full-rank-preserving swap hill-climb (accept non-worsening fitness)
        cur_fit = _fit(d0, fr0, ad0)
        cur_set = set(cur_pts)
        for _step in range(climb_steps):
            outside = [x for x in allpts if x not in cur_set]
            rng.shuffle(outside)
            drop_idx = rng.randrange(len(cur_pts))
            base = cur_pts[:drop_idx] + cur_pts[drop_idx + 1:]  # still a cap (subset)
            moved = False
            for newx in outside[:swap_tries]:
                if not cap_extends(base, newx, p):
                    continue
                cand_pts = base + [newx]
                dd, ffr, ad_new, _ = _eval(points_to_columns(cand_pts, p))
                cand_fit = _fit(dd, ffr, ad_new)
                if cand_fit >= cur_fit:
                    cur_pts, cur_set, cur_fit = cand_pts, set(cand_pts), cand_fit
                    moved = True
                    break
            if not moved:
                break
    return best


def random_search(p, m, trials, seed=0, target_k=None, max_distance=6, n_jobs=1,
                  sampler="uniform", climb_steps=30, swap_tries=8) -> list:
    """Randomized puncture-location search over RM_p(r_max,m) (NOTES sec 5).

    Evaluates ``trials`` candidate puncture sets, keeping every full-rank, distance-certified
    code (deduped by puncture set), sorted by gamma. ``max_distance`` bounds the (certified)
    distance; codes above it are skipped, not mis-reported. ``target_k`` fixes the puncture
    count per candidate.

    ``sampler`` (see SAMPLERS): "uniform" (default, unchanged i.i.d. sampling); "capset"
    (draw cap sets -- no 3 collinear points -- reproduces the cap-structured paper codes uniform
    misses); "plane_spread" (draw caps that ALSO have no 4 coplanar points, i.e. <=3 per 2-flat --
    reaches the higher-distance codes the plain cap misses, e.g. [[230,13,6]]_3 d=6 vs the cap's
    d=5); "capset_climb" (cap seed + cap-preserving distance hill-climb); "arc_climb" (cap seed +
    a hill-climb on the lexicographic (distance, -A_d) fitness -- targets d>=4 by driving the
    minimum-weight-codeword multiplicity to zero where plain distance is a flat plateau; needs
    the exact A_d engine, so it only gets its gradient in the small-dual regime, otherwise it
    degrades to capset_climb). The structure-aware samplers need ``target_k`` within the
    (plane-spread) cap-size bound to be useful; ``climb_steps``/``swap_tries`` tune the climb.
    "capset"/"plane_spread" are seed-only (one draw per trial); the climb samplers do a full
    seed+climb per trial.

    ``n_jobs`` controls process-level parallelism (trials are independent): n_jobs=1 is serial
    and deterministic from ``seed``; n_jobs>1/-1 splits trials across worker processes (joblib),
    each with a seed-derived RNG. Parallel results are reproducible for a fixed
    (seed, sampler, n_jobs, trials) but are not bit-identical to the serial stream.
    """
    if sampler not in SAMPLERS:
        raise ValueError(f"sampler must be one of {SAMPLERS}, got {sampler!r}")
    pm = p ** m
    r = r_max(m, p)
    cap = min(pm - 1, max(2, 2 * d_rm(r, m, p)))
    args = (p, m, r, cap, pm)
    tail = (target_k, max_distance, sampler, climb_steps, swap_tries)

    if n_jobs == 1:
        merged = _search_chunk(*args, trials, seed, *tail)
    else:
        from joblib import Parallel, delayed
        n = (os.cpu_count() or 1) if n_jobs in (-1, None) else n_jobs
        n = max(1, min(n, trials))
        per = [trials // n + (1 if i < trials % n else 0) for i in range(n)]
        chunks = Parallel(n_jobs=n)(
            delayed(_search_chunk)(*args, per[i], seed * 100003 + i, *tail)
            for i in range(n) if per[i] > 0
        )
        merged = {}
        for d in chunks:
            merged.update(d)

    out = list(merged.values())
    out.sort(key=lambda c: c.gamma)
    return out


def _cost(c: Code, delta_in: float) -> float:
    """Single-round distillation cost C = n / nbar_T (NOTES eq 39)."""
    return cost(c.n, nbar_T(c.n, c.k, c.p, delta_in))


def search(p, m_values=None, trials_per_m=2000, seed=0, delta_in=1e-3, top=10, n_jobs=1,
           sampler="uniform") -> dict:
    """Top-level: given a prime p, redo the paper's search across a range of m.

    Runs the analytic Manhattan sweep (all m) plus a randomized explicit search (small p**m),
    and returns the best codes ranked by yield gamma and by single-round cost C. ``n_jobs`` and
    ``sampler`` are passed to the explicit ``random_search``. Returns {best_by_gamma,
    best_by_cost, scanned}.
    """
    assert isinstance(p, int) and p >= 2
    if m_values is None:
        m_values = [m for m in range(2, 12) if p ** m <= 2000] or [2, 3]

    found: list[Code] = []
    for m in m_values:
        found.extend(manhattan_sweep(p, m))
        if p ** m <= EXPLICIT_MAX_BLOCK:
            found.extend(random_search(p, m, trials_per_m, seed=seed, n_jobs=n_jobs,
                                       sampler=sampler))

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
