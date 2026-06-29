"""Structured *punctured minimum distance* (d, not full A_d) of m=2 punctured RM codes.

For a punctured triorthogonal code the quantum code ``G0^perp`` is the *punctured*
generalized Reed-Muller code ``RM_p(rtilde, m)`` restricted to the surviving
coordinates (``rtilde = m(p-1) - r - 1``).  When the punctured-column submatrix is full
rank the restriction map ``RM_p(rtilde,m) -> G0^perp`` is a bijection, so

    d = min_{c != 0 in RM_p(rtilde,m)} |supp(c) \\ S|        (S = puncture set, |S| = k).

This module pins ``d`` for ``m = 2`` *cheaply*, sidestepping the ``(p-1)^a * C(n,a)``
meet-in-the-middle blow-up that makes the direct weight-d search ~9e9 ops at ``p = 17``.

The geometric structure (m = 2, alpha = 1)
------------------------------------------
Write ``rtilde = alpha (p-1) + beta``.  For ``alpha = 1`` (the regime of the new p=17,
p=19 search codes, ``(p-1) <= rtilde < 2(p-1)``) the **minimum-weight codewords of
``RM_p(rtilde,2)`` are supported on a single affine line** of ``AG(2,p)`` (Delsarte-
Goethals-MacWilliams): on a line ``ell`` the codewords supported in ``ell`` are exactly the
univariate polynomials of degree ``<= deg_line := rtilde - (p-1) = beta`` (the restricted
code ``RM_p(beta, 1)``, an MDS ``[p, beta+1, p-beta]`` Reed-Solomon code).  A nonzero
degree-``<=beta`` polynomial has at most ``beta`` roots, so on a line ``ell`` with
``s = |S cap ell|`` punctured points (``T = ell \\ S`` survivors) the **minimum number of
surviving non-zeros is exactly ``|T| - beta = (p - s) - beta``** (place all ``beta`` roots
on survivors).  Hence the line-supported punctured distance is

    d_lines = (p - beta) - max_ell |S cap ell|          [ = d_RM - max_ell |S cap ell| ].

There are only ``p(p+1)`` lines, so ``max_ell |S cap ell|`` is a trivial scan -- no
``C(n,d)`` enumeration.  ``d_lines`` is a **certified upper bound** on ``d`` (we emit an
explicit weight-``d_lines`` codeword and verify ``X_stab . v = 0``), and it equals the true
``d`` whenever no genuinely 2-dimensional codeword undercuts it -- which is confirmed by the
certified MITM lower bound ``qmsd.mindist.min_dependent_columns`` ruling out every codeword
(line *or* full-span) of weight ``< d_lines``.  Critically, that lower-bound search only has
to reach ``d_lines - 1`` (e.g. weight <=5 for d=6), so it AVOIDS the expensive weight-6
right-half MITM that exhausts >10 min -- the structured upper bound supplies the d<=6 half.

Validated: reproduces ``d = 5`` for ``[[234,55,5]]_17`` and pins ``d = 6`` for
``[[237,52,>=6]]_17`` with explicit codewords (see ``tests/test_structured_pe.py``).
"""
from __future__ import annotations

import numpy as np

from .reedmuller import r_tilde, d_rm
from .triorthogonal import build_triorthogonal_code

__all__ = ["line_punctured_distance", "enumerate_lines"]


def enumerate_lines(p):
    """Yield ``(point_indices, t_params)`` for each of the ``p(p+1)`` affine lines of AG(2,p).

    Point index convention matches ``qmsd.reedmuller.points``: column ``c`` <-> point
    ``(x, y) = (c % p, c // p)`` (x least significant), 0-indexed.  Each line is returned as
    its ``p`` global point indices together with the natural line parameter ``t`` of each
    point (``t = y`` for the vertical lines ``x = a``; ``t = x`` for the sloped lines
    ``y = s x + b``), so a univariate codeword in ``t`` lifts directly.
    """
    rng = list(range(p))
    # vertical lines x = a, parametrised by t = y
    for a in range(p):
        yield [a + p * t for t in rng], rng
    # sloped lines y = s*x + b, parametrised by t = x
    for s in range(p):
        for b in range(p):
            yield [x + p * ((s * x + b) % p) for x in rng], rng


def _line_codeword(line_pts, t_params, S, surv_index, beta, p, n):
    """Explicit weight-(|T|-beta) codeword of ``G0^perp`` supported on a single line.

    Builds the univariate degree-``beta`` polynomial ``g(t) = prod (t - root_i)`` whose
    ``beta`` roots are the first ``beta`` surviving points of the line, evaluates it on the
    line's surviving points, and zeros it everywhere else.  By the structured identity this
    is the restriction to the surviving coordinates of an ``RM_p(rtilde,2)`` codeword, i.e.
    a genuine codeword of ``G0^perp = ker(X_stab)``.  Returns a length-``n`` int vector.
    """
    surv_on = [(pt, t) for pt, t in zip(line_pts, t_params) if pt not in S]
    if len(surv_on) < beta:
        raise ValueError("line has fewer survivors than beta; cannot place all roots")
    roots = [t for _pt, t in surv_on[:beta]]

    def g(t):
        v = 1
        for rt in roots:
            v = (v * ((t - rt) % p)) % p
        return v

    v = np.zeros(n, dtype=np.int64)
    for pt, t in zip(line_pts, t_params):
        if pt not in S:
            v[surv_index[pt]] = g(t) % p
    return v


def line_punctured_distance(
    p,
    m,
    r,
    puncture_columns_1indexed,
    lower_bound=None,
    certify_lower=False,
    G=None,
):
    """Structured (line-based) punctured minimum distance of the m=2 punctured-RM code.

    Parameters
    ----------
    p, m, r : int
        Prime ``p``, number of RM variables ``m`` (must be 2), base degree ``r`` (= r_max).
        ``rtilde = m(p-1) - r - 1`` must satisfy ``p-1 <= rtilde < 2(p-1)`` (``alpha = 1``),
        the regime where minimum-weight ``RM_p(rtilde,2)`` codewords are line-supported.
    puncture_columns_1indexed : iterable[int]
        Puncture set ``S`` (1-indexed; column ``c`` <-> point ``c-1``).
    lower_bound : int or None
        A known certified lower bound on the true distance (e.g. from a prior
        ``min_dependent_columns(X_stab, p, d_max=d_lines-1)`` run).  If ``>= d_lines`` the
        result is marked exact.
    certify_lower : bool
        If True, run ``min_dependent_columns(X_stab, p, d_max=d_lines-1)`` to certify the
        lower bound here.  WARNING: at ``p = 17`` the weight-``(d_lines-1)`` MITM can take
        many minutes; default False -- pass a precomputed ``lower_bound`` instead.

    Returns
    -------
    dict with keys
        ``d_lines``            : line-supported punctured distance (certified UPPER bound on d).
        ``distance``           : the exact d when a lower bound confirms it, else None.
        ``exact``              : bool -- whether ``distance`` is pinned.
        ``d_RM``               : minimum distance of ``RM_p(rtilde,2)`` (= p - beta).
        ``beta``               : ``rtilde - (p-1)`` (univariate restricted degree on a line).
        ``max_line_punctures`` : ``max_ell |S cap ell|``.
        ``best_line_points``   : the achieving line's global point indices (0-indexed).
        ``certificate``        : {weight, support, values, in_kernel} of the explicit codeword.
        ``lower_bound``        : the lower bound used (passed or certified), if any.
    """
    assert isinstance(p, int) and p >= 2
    if m != 2:
        raise NotImplementedError(
            "line_punctured_distance handles m=2 only (minimum-weight RM_p(rtilde,2) "
            "codewords are line-supported); use qmsd.structured_ad for the general flat "
            "decomposition."
        )
    rtilde = r_tilde(r, m, p)
    beta = rtilde - (p - 1)
    if not (0 <= beta <= p - 1):
        raise NotImplementedError(
            f"rtilde={rtilde} is outside the alpha=1 regime (need {p-1} <= rtilde < "
            f"{2*(p-1)}); line-supported minimum-weight structure does not apply."
        )
    d_RM = d_rm(rtilde, m, p)
    assert d_RM == p - beta, "internal: d_RM mismatch for alpha=1"

    code = build_triorthogonal_code(p, m, r, puncture_columns_1indexed, G=G)
    if not code["full_rank"]:
        raise ValueError("puncture submatrix not full rank; RM<->G0^perp bijection fails")
    X = np.asarray(code["X_stab"]).astype(np.int64) % p
    n = code["n"]

    S = set(int(c) - 1 for c in puncture_columns_1indexed)
    surv = [c for c in range(p ** m) if c not in S]
    assert len(surv) == n
    surv_index = {c: i for i, c in enumerate(surv)}

    # --- the only loop: count punctures on each of the p(p+1) lines ---
    max_pun = -1
    best_line = None
    for line_pts, t_params in enumerate_lines(p):
        s = sum(1 for pt in line_pts if pt in S)
        if s > max_pun:
            max_pun = s
            best_line = (line_pts, t_params)

    d_lines = d_RM - max_pun  # = (p - beta) - max_ell |S cap ell|
    if d_lines < 1:
        raise ValueError(
            f"d_lines = {d_lines} < 1: a line carries a (near-)zero-weight punctured "
            f"codeword (max |S cap line| = {max_pun} >= d_RM = {d_RM}); the code is degenerate."
        )

    # --- explicit certificate codeword of weight d_lines on the heaviest line ---
    line_pts, t_params = best_line
    v = _line_codeword(line_pts, t_params, S, surv_index, beta, p, n)
    weight = int(np.count_nonzero(v))
    in_kernel = bool(np.all((X @ v) % p == 0))
    support = [c for c in surv if v[surv_index[c]] != 0]
    values = {int(c): int(v[surv_index[c]]) for c in support}
    certificate = {
        "weight": weight,
        "support": support,            # 0-indexed surviving point columns
        "values": values,
        "in_kernel": in_kernel,        # X_stab . v == 0  (genuine logical operator)
    }
    # the certificate must be a genuine weight-d_lines logical operator
    assert in_kernel, "explicit line codeword is not in ker(X_stab) -- structure violated"
    assert weight == d_lines, f"explicit codeword weight {weight} != d_lines {d_lines}"

    # --- lower bound: certify d >= d_lines (rules out line *and* 2D codewords < d_lines) ---
    lb = lower_bound
    if certify_lower:
        from .mindist import min_dependent_columns
        try:
            got = min_dependent_columns(X, p, d_max=d_lines - 1)
            lb = got  # a smaller dependency exists: lines were NOT tight, true d = got
        except ValueError:
            lb = max(lb or 0, d_lines)  # no dependency <= d_lines-1  =>  d >= d_lines

    exact = lb is not None and lb >= d_lines
    distance = d_lines if exact else None

    return {
        "d_lines": d_lines,
        "distance": distance,
        "exact": exact,
        "d_RM": d_RM,
        "beta": beta,
        "max_line_punctures": max_pun,
        "best_line_points": list(line_pts),
        "certificate": certificate,
        "lower_bound": lb,
    }
