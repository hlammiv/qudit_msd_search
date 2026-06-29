"""Code dataclass and constructors for the puncture / Manhattan families (NOTES sec 5, sec 7)."""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .reedmuller import r_max, r_tilde, d_rm
from .pnomial import pnomial_le, pnomial_gt
from .distance import delta_p, A_d_logical_Z
from .triorthogonal import build_triorthogonal_code
from .mindist import min_dependent_columns


@dataclass(frozen=True)
class Code:
    """A [[n,k,d]]_p triorthogonal distillation code record (NOTES sec 7)."""

    p: int
    n: int
    k: int
    d: int | None
    m: int | None = None
    r: int | None = None
    w: int | None = None
    puncture_columns: tuple | None = None
    A_d: int | None = None
    full_rank: bool | None = None
    d_certified: bool = True
    A_d_exact: bool | None = None  # True when A_d came from the exact MacWilliams engine

    @property
    def gamma(self) -> float:
        """gamma = log(n/k)/log(d), the distillation overhead exponent (NOTES sec 7).

        Returns NaN when the distance is unknown/uncertified or too small to
        define a yield (d is None or d <= 1, or k <= 0).
        """
        if self.d is None or self.d <= 1 or self.k <= 0 or self.n <= 0:
            return float("nan")
        return math.log(self.n / self.k) / math.log(self.d)

    @property
    def label(self) -> str:
        d = "?" if self.d is None else self.d
        return f"[[{self.n},{self.k},{d}]]_{self.p}"


def gamma(n, k, d) -> float:
    """gamma = log(n/k)/log(d) (NOTES sec 7)."""
    return math.log(n / k) / math.log(d)


# Default enumeration budget (number of messages q**dim(G0)) for the exact MacWilliams
# path. ``code_from_puncture`` keeps it OFF by default (exact_budget=0) so the existing
# behavior -- and the performance-sensitive search.py -- are byte-for-byte unchanged.
# ``code_certify`` turns it on with this budget.
DEFAULT_EXACT_BUDGET = 5_000_000


def code_from_puncture(
    p, m, puncture_columns_1indexed, r=None, compute_A_d=True, max_distance=None, G=None,
    exact_budget=0,
) -> Code:
    """Build a Code from an explicit puncture-column set (NOTES sec 5, sec 7).

    r defaults to r_max(m,p). The triorthogonal space RM_p(r,m) is punctured at the
    1-indexed columns; the quantum distance is d(G0^perp) = min weight of the dual of
    the X-stabilizer G0 (NOTES Thm 3). A_d is computed when feasible.

    If ``max_distance`` is given and the true distance exceeds it, the distance is left
    UNCERTIFIED (d=None, d_certified=False) rather than under-reported -- callers (e.g.
    the search) should skip such codes. ``max_distance`` is also what keeps the distance
    search tractable for larger blocks.

    ``exact_budget`` (default 0 = OFF) enables the EXACT MacWilliams engine
    (``weightdist.exact_distance_and_Ad``) when ``p**dim(G0) <= exact_budget``: it
    enumerates the small generator ``G0`` and MacWilliams-transforms to the FULL weight
    distribution of ``G0^perp``, certifying the distance for ANY d (lifting the d<=6 MITM
    cap) and the exact ``A_d`` (closing the large-code A_d gap), in the small-dual /
    high-puncture regime. When the exact path is infeasible or declines, it falls back to
    the meet-in-the-middle ``min_dependent_columns`` / ``A_d_logical_Z`` exactly as before.
    Default 0 keeps every existing caller and test unchanged.
    """
    if r is None:
        r = r_max(m, p)
    built = build_triorthogonal_code(p, m, r, puncture_columns_1indexed, G=G)
    n, k, full_rank = built["n"], built["k"], built["full_rank"]

    G0 = np.asarray(built["X_stab"]) % p
    d: int | None = None
    d_certified = True
    A_d = None
    A_d_exact: bool | None = None

    # --- Exact MacWilliams path (opt-in via exact_budget) ----------------------------
    # When dim(G0) is small enough to enumerate, the dual weight distribution is exact, so
    # both the distance (uncapped) and A_d are certified. The minimum-weight dual codewords
    # are all LOGICAL here (no weight-d stabilizers; empirically certified for every Table-3
    # code, see weightdist + tests), so A_d_logical == B_d.
    if exact_budget:
        from .weightdist import exact_distance_and_Ad
        res = exact_distance_and_Ad(G0, p, max_words=exact_budget)
        if res["feasible"] and res["distance"] is not None:
            d = res["distance"]
            d_certified = True
            if compute_A_d:
                A_d = res["A_d_logical"]
                A_d_exact = True

    # --- Meet-in-the-middle fallback (capped at d<=6) --------------------------------
    if d is None:
        # Quantum distance = min weight of G0^perp = min number of F_p-dependent columns of
        # G0 (= X_stab, the parity check of G0^perp), via the meet-in-the-middle column
        # search (NOTES Thm 3). Certified exactly for distance up to 6.
        try:
            d = min_dependent_columns(G0, p, d_max=max_distance)
        except (ValueError, OverflowError):
            # distance exceeds the cap / not certifiable within budget
            d, d_certified = None, False

        if compute_A_d and d is not None:
            try:
                A_d = A_d_logical_Z(built, p, d)
                A_d_exact = True
            except NotImplementedError:
                A_d = None  # exact A_d infeasible for this block size; left unknown

    return Code(
        p=p, n=n, k=k, d=d, m=m, r=r, w=None,
        puncture_columns=tuple(puncture_columns_1indexed),
        A_d=A_d, full_rank=full_rank, d_certified=d_certified, A_d_exact=A_d_exact,
    )


def code_certify(
    p, m, puncture_columns_1indexed, r=None, compute_A_d=True, G=None,
    exact_budget=DEFAULT_EXACT_BUDGET,
) -> Code:
    """Build a Code certifying distance + A_d via the EXACT MacWilliams engine when feasible.

    Thin wrapper over ``code_from_puncture`` with the exact path enabled by default
    (``exact_budget=DEFAULT_EXACT_BUDGET``). In the small-dual / high-puncture regime
    (``p**dim(G0) <= exact_budget``) this certifies the distance for ANY ``d`` -- lifting the
    ``d<=6`` cap of the meet-in-the-middle search -- and the exact ``A_d`` even for large
    blocks where ``A_d_logical_Z`` would refuse (its ``C(n,d)`` subset scan blows the budget).
    Outside that regime it falls back to the meet-in-the-middle distance / ``A_d_logical_Z``.

    ``Code.A_d_exact`` reports whether ``A_d`` is an exact certified count.
    """
    return code_from_puncture(
        p, m, puncture_columns_1indexed, r=r, compute_A_d=compute_A_d, max_distance=None,
        G=G, exact_budget=exact_budget,
    )


def code_from_manhattan(p, m, w, r=None) -> Code:
    """Build the analytic Manhattan-family Code, exact via Theorem 4 (NOTES Thm 4, 5; sec 5).

    No matrices: n = [m,>w]_p, k = [m,<=w]_p, and the distance is the closed-form
    Delta_p(m, rtilde, w) evaluated at the dual degree rtilde = m(p-1)-r-1. This scales
    to the astronomically large Table-2 codes. r defaults to r_max(m,p).
    """
    if r is None:
        r = r_max(m, p)
    rt = r_tilde(r, m, p)
    n = pnomial_gt(m, w, p)
    k = pnomial_le(m, w, p)
    d = delta_p(m, rt, w, p)
    # Full rank is guaranteed when the number of punctures stays below d_RM (NOTES eq 15 fn 3).
    full_rank = (k < d_rm(r, m, p)) or None
    return Code(
        p=p, n=n, k=k, d=d, m=m, r=r, w=w,
        puncture_columns=None, A_d=None, full_rank=full_rank, d_certified=True,
    )
