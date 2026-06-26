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


def code_from_puncture(
    p, m, puncture_columns_1indexed, r=None, compute_A_d=True, max_distance=None, G=None
) -> Code:
    """Build a Code from an explicit puncture-column set (NOTES sec 5, sec 7).

    r defaults to r_max(m,p). The triorthogonal space RM_p(r,m) is punctured at the
    1-indexed columns; the quantum distance is d(G0^perp) = min weight of the dual of
    the X-stabilizer G0 (NOTES Thm 3). A_d is computed when feasible.

    If ``max_distance`` is given and the true distance exceeds it, the distance is left
    UNCERTIFIED (d=None, d_certified=False) rather than under-reported -- callers (e.g.
    the search) should skip such codes. ``max_distance`` is also what keeps the distance
    search tractable for larger blocks.
    """
    if r is None:
        r = r_max(m, p)
    built = build_triorthogonal_code(p, m, r, puncture_columns_1indexed, G=G)
    n, k, full_rank = built["n"], built["k"], built["full_rank"]

    # Quantum distance = min weight of G0^perp = min number of F_p-dependent columns of
    # G0 (= X_stab, the parity check of G0^perp), via the meet-in-the-middle column search
    # (NOTES Thm 3). Certified exactly for distance up to 6 (covers every Table-3 code).
    G0 = np.asarray(built["X_stab"]) % p
    d: int | None
    d_certified = True
    try:
        d = min_dependent_columns(G0, p, d_max=max_distance)
    except (ValueError, OverflowError):
        # distance exceeds the cap / not certifiable within budget
        d, d_certified = None, False

    A_d = None
    if compute_A_d and d is not None:
        try:
            A_d = A_d_logical_Z(built, p, d)
        except NotImplementedError:
            A_d = None  # exact A_d infeasible for this block size; left unknown

    return Code(
        p=p, n=n, k=k, d=d, m=m, r=r, w=None,
        puncture_columns=tuple(puncture_columns_1indexed),
        A_d=A_d, full_rank=full_rank, d_certified=d_certified,
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
