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
empirically (max-2flat 3 -> d=6, 4 -> d=5, ...). Phase 2 (the certified d>=7 LOWER bound) needs
the weight-class enumeration up to weight w+k and inherits the full-span crux; not in scope here.
"""
from __future__ import annotations

from itertools import product

import numpy as np

from .reedmuller import rm_generator, d_rm
from .structured_ad import _flats

__all__ = ["structured_distance", "max_flat_occupancy"]


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
