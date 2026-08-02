"""Geometric (flat-occupancy) distance UPPER bound / fast screen for punctured-RM codes.

Distance = min over dual RM codewords of |supp(c) \\ S| (S = puncture set). RM min-weight
codewords are supported on affine flats (Delsarte-Goethals-MacWilliams / Leducq), so
enumerating flat-supported codewords of affine span <= ``jmax`` and taking the minimum
punctured weight gives an UPPER bound on the distance:

    d(S) <= d_geo = min over flat codewords c of |supp(c) \\ S|.

It is EXACT whenever the distance-binding codeword is flat-supported (the common case for the
paper's high-distance arc codes; the residual is the full-span crux, D_CRUX_REDUCTION.md).
There is no C(n,d) column scan and no p^(d/2) weight table, so it reaches the d>=7 regime the
meet-in-the-middle certifier (capped at d<=6) and Brouwer-Zimmermann cannot.

Phase 1 (this module): the upper bound / screen. For the paper's qutrit m=5 codes the j=2
term is exactly  d_RM - max_{2-flat} |plane cap S|  -- the flat-occupancy law observed
empirically (max-2flat 3 -> d=6, 4 -> d=5, ...).

Phase 2 (weight_hierarchy / flat_lower_bound): the certified LOWER bound. Its headline result is
a rigorous NO-GO, not a d>=7 certificate: at qutrit m=5, ANY 3 puncture points lie on a common
2-flat (they span <= 2 dims), whose weight-9 indicator is a codeword, so |supp\\S| <= 9 - 3 = 6.
Hence  d(S) <= 6 for EVERY puncture set of size k>=3  -- so [[230,13,6]] is distance-OPTIMAL and
d>=7 is IMPOSSIBLE at m=5. In the small-k regime the lower bound meets the Phase-1 upper bound and
pins the distance EXACTLY with no MITM. The general large-k lower bound still inherits the
full-span crux (D_CRUX_REDUCTION.md).
"""
from __future__ import annotations

from itertools import product

import numpy as np

from .reedmuller import rm_generator, d_rm
from .structured_ad import _flats

__all__ = ["structured_distance", "max_flat_occupancy", "geometric_distance_upper",
           "geometric_distance_dual", "weight_hierarchy", "flat_lower_bound"]


def _min_punctured_weight(Gj, surv_mask, p, dim_cap=12):
    """Min over nonzero codewords c = msg @ Gj of the number of SURVIVING nonzeros of c.

    Gj is the (dim x p^j) generator of the flat-restricted code C_F. For the min-weight term
    dim is tiny (j=2, r'=0 -> dim 1, the flat indicator). Returns None if dim exceeds dim_cap
    (that flat's codewords are higher-weight and do not bind the minimum -- skipped)."""
    G = np.asarray(Gj, dtype=np.int64) % p
    dim = G.shape[0]
    if dim == 0 or dim > dim_cap:
        return None
    best = None
    for msg in product(range(p), repeat=dim):
        if not any(msg):
            continue
        c = (np.asarray(msg, dtype=np.int64) @ G) % p
        w = int(np.count_nonzero(c[surv_mask]))
        if w and (best is None or w < best):
            best = w
    return best


def structured_distance(p, m, r, puncture_columns, jmax=2):
    """Geometric UPPER bound on the punctured distance from flat-supported codewords (span<=jmax).

    Returns a dict:
        ``d_upper``   : min punctured weight over the enumerated flat codewords (<= true d).
        ``d_RM``      : min weight of the capping code RM_p(rtilde,m).
        ``exact_if``  : note on when d_upper is the true distance.
        ``jmax``      : cap used.
    Exact when the distance-binding codeword is flat-supported (jmax=m, or empirically for the
    high-distance arc codes). This is a certified upper bound regardless (an explicit codeword)."""
    rtilde = m * (p - 1) - r - 1
    S = set(int(c) - 1 for c in puncture_columns)
    d_min_rm = d_rm(rtilde, m, p)
    best = None
    for j in range(2, min(jmax, m) + 1):
        rr = rtilde - (m - j) * (p - 1)
        if rr < 0:
            continue
        Gj = np.asarray(rm_generator(rr, j, p), dtype=np.int64) % p   # (dim, p^j), col t <-> flat point t
        for colmap in _flats(m, j, p):
            surv_mask = np.fromiter((int(c) not in S for c in colmap), dtype=bool, count=len(colmap))
            if not surv_mask.any():
                continue
            w = _min_punctured_weight(Gj, surv_mask, p)
            if w is not None and (best is None or w < best):
                best = w
    return {
        "d_upper": best,
        "d_RM": d_min_rm,
        "jmax": min(jmax, m),
        "exact_if": "binding codeword is flat-supported (span<=jmax); true d otherwise strictly less",
    }


def geometric_distance_dual(p, m, r, puncture_columns, j=2, d_max=6):
    """Point-restricted geometric distance UPPER bound via the flat-restricted DUAL MITM -- works
    where structured_distance's codeword-enumeration OOMs (high-dim flats, e.g. p=7 m=4 where the
    j=2 restricted code RM_7(4,2) has dim 15 >> the enum cap).

    For each affine j-flat F, the flat-restricted dual code is RM_p(rr, j), rr = rtilde-(m-j)(p-1);
    its min punctured weight (= min over dual codewords supported on F of |supp\\S|) is the min
    distance of that code restricted to the surviving (F\\S) columns, computed by
    min_dependent_columns(dual_matrix(G[:, F\\S]), p). Only flats with |S cap F| >= d_RM - d_max can
    bind at <= d_max (min punctured >= d_RM - |S cap F|), so lighter flats are skipped. Returns the
    min over flats (or None if none bind within d_max). Validated EXACT on the qutrit m=5 oracle
    codes; at p=7 m=4 it returns None on the random d=5 codes -> binding codeword NOT j=2-supported
    (full-span crux). This is a certified UPPER bound; == true d iff the binding codeword is
    flat-supported (span <= j)."""
    from .reedmuller import rm_generator
    from .triorthogonal import dual_matrix
    from .mindist import min_dependent_columns
    from .structured_ad import _flats
    rtilde = m * (p - 1) - r - 1
    rr = rtilde - (m - j) * (p - 1)
    if rr < 0:
        return None
    G = np.asarray(rm_generator(rr, j, p), dtype=np.int64) % p
    d_min_rm = d_rm(rtilde, m, p)
    S = set(int(c) - 1 for c in puncture_columns)
    min_occ = d_min_rm - d_max
    best = None
    for colmap in _flats(m, j, p):
        surv = np.fromiter((int(c) not in S for c in colmap), dtype=bool, count=len(colmap))
        if int((~surv).sum()) < min_occ:          # |S cap F| too small to bind within d_max
            continue
        surv_cols = np.where(surv)[0]
        try:
            w = min_dependent_columns(dual_matrix(G[:, surv_cols], p), p, d_max=d_max)
        except ValueError:
            continue                              # min punctured weight > d_max on this flat
        if best is None or w < best:
            best = w
    return best


def geometric_distance_upper(p, m, r, puncture_columns):
    """Fast CERTIFIED UPPER bound on the distance from the full-2-flat codeword, or ``None`` when
    that codeword does not exist for this (p,m,r) (caller then must NOT screen).

    A full affine 2-flat carries p^2 points; its indicator is a product of (m-2) linear-form
    factors, degree (m-2)(p-1), so it lies in the dual RM_p(rtilde,m) iff (m-2)(p-1) <= rtilde.
    When it does, its punctured weight is p^2 - |flat cap S|, so

        d(S) <= p^2 - max_2flat_occupancy(S)

    -- one flat-incidence pass (~tens of ms), no MITM. Used as a distance FLOOR pre-screen: if
    this bound is below the floor, the true distance is too, so the ~15s weight-3 MITM is skipped.
    It is TIGHT (== true d) exactly when the 2-flat is the minimum-weight support, i.e. rtilde ==
    (m-2)(p-1) (p^2 == d_RM) -- the qutrit m=5 regime (9 - max-2flat). Elsewhere it is a valid but
    LOOSER bound (rarely below the floor -> few skips, but never an unsound drop). Returns None
    when (m-2)(p-1) > rtilde so the caller falls back to the MITM rather than trust a non-codeword."""
    rtilde = m * (p - 1) - r - 1
    if (m - 2) * (p - 1) > rtilde:
        return None
    return p * p - max_flat_occupancy(p, m, puncture_columns, 2)


def max_flat_occupancy(p, m, puncture_columns, j=2):
    """Max number of puncture points on any affine j-flat -- the geometric quantity that, for
    the qutrit m=5 codes, gives d = d_RM - max_2flat_occupancy. Pure incidence geometry."""
    S = set(int(c) - 1 for c in puncture_columns)
    best = 0
    for colmap in _flats(m, j, p):
        occ = sum(1 for c in colmap if int(c) in S)
        if occ > best:
            best = occ
    return best


# --- Phase 2: certified LOWER bound (weight-gap) ------------------------------------------
# Distance splits over the RM weight classes. The min-weight (span-2) codewords give an
# EXHAUSTIVE term d_RM - max_2flat_occupancy(S). Every heavier codeword has RM-weight >= w2
# (the second weight), so its punctured weight is >= w2 - k. Hence
#       d(S) >= min( d_RM - max_2flat_occupancy(S),  w2 - k ).
# NO-GO COROLLARY (m=5): any 3 points share a 2-flat => max_2flat_occupancy(S) >= 3 for every
# k>=3 => d_span2 = d_RM - max_2flat <= 9 - 3 = 6. So d(S) <= 6 for ALL k>=3: d>=7 is impossible
# and [[230,13,6]] is distance-optimal at qutrit m=5. (This is the Phase-1 upper bound applied
# structurally; the lower bound below pins the exact value when it meets that upper bound.)
# This is a certified lower bound; combined with the Phase-1 upper bound it pins the distance
# EXACTLY when they meet -- reaching d>=7 with no MITM in the small-k / short-window regime.
# SOUND iff w2 is a true lower bound on every codeword weight above d_RM. For RM_3(6,5) the
# minimal-affine-span hierarchy gives span-2=9, span-3=12; confirming w2=12 for ALL spans
# (no weight-10/11 codeword of span 4,5) is the Leducq weight-hierarchy question = the
# full-span crux (see GEOMETRIC_CERTIFIER_SCOPE.md). Pass w2 explicitly.
def _fullspan_min_weight(p, j, rr, dim_cap=13):
    """Min weight of a codeword of RM_p(rr,j) whose affine span is exactly j (brute; small dim)."""
    G = np.asarray(rm_generator(rr, j, p), dtype=np.int64) % p
    dim = G.shape[0]
    if dim > dim_cap:
        return None
    pts = np.array(list(product(range(p), repeat=j)), dtype=np.int64)
    best = None
    for msg in product(range(p), repeat=dim):
        if not any(msg):
            continue
        c = (np.asarray(msg, dtype=np.int64) @ G) % p
        supp = np.nonzero(c)[0]
        w = supp.size
        if best is not None and w >= best:
            continue
        P = (pts[supp] - pts[supp][0]) % p
        if _fp_rank(P, p) == j:            # spans the full j-flat
            best = w
    return best


def _fp_rank(M, p):
    M = np.asarray(M, dtype=np.int64) % p
    nr, nc = M.shape
    rk = 0
    for c in range(nc):
        piv = next((i for i in range(rk, nr) if M[i, c] % p), None)
        if piv is None:
            continue
        M[[rk, piv]] = M[[piv, rk]]
        M[rk] = (M[rk] * pow(int(M[rk, c]), p - 2, p)) % p
        for i in range(nr):
            if i != rk and M[i, c] % p:
                M[i] = (M[i] - M[i, c] * M[rk]) % p
        rk += 1
        if rk == nr:
            break
    return rk


def weight_hierarchy(p, m, r, jmax=3):
    """{span j : min weight of an exactly-span-j codeword} of RM_p(rtilde,m), for j=2..jmax
    (brute per flat; small spans only). The second weight w2 is min over j>=3."""
    rtilde = m * (p - 1) - r - 1
    out = {}
    for j in range(2, min(jmax, m) + 1):
        rr = rtilde - (m - j) * (p - 1)
        if rr < 0:
            continue
        w = _fullspan_min_weight(p, j, rr)
        if w is not None:
            out[j] = w
    return out


def flat_lower_bound(p, m, r, puncture_columns, w2):
    """Certified lower bound d(S) >= min( d_RM - max_2flat_occupancy(S), w2 - k ), with w2 a
    (sound) lower bound on every RM-weight above d_RM. Meets the Phase-1 upper bound -> EXACT
    distance, no MITM, in the small-k regime (incl. d>=7). See module notes for the w2 caveat."""
    S = set(int(c) - 1 for c in puncture_columns)
    k = len(S)
    d_RM = d_rm(m * (p - 1) - r - 1, m, p)
    d_span2 = d_RM - max_flat_occupancy(p, m, puncture_columns, 2)
    return min(d_span2, w2 - k)
