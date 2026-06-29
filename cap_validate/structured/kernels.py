"""Numba hot loops for Method 2.

Two kernels:

1. ``max_topj`` -- the heavy one.  Given the syndrome (coset) maps of all directions
   of a fixed dimension and the cap points, compute, per direction, the top-j
   coset-occupancy sum  max_{j distinct cosets} |union cap S|  and return the whole
   per-direction array (so the caller takes the global max -> rigorous max|supp cap S|
   for the corresponding "union of j parallel flats" family).  Parallel over directions.

2. ``inter_counts`` -- a brute |supp cap S| over explicit support index-sets, used only
   as an independent cross-check of the additive ``max_topj`` result on a sample.
"""
from __future__ import annotations

import numpy as np
from numba import njit, prange

P = 3


@njit(cache=True, parallel=True, fastmath=False)
def max_topj(Hs, Spts, codim, j):
    """Per-direction top-j coset-occupancy sum.

    Hs   : (Nsub, codim, 7) int64 syndrome maps (ker H = direction V).
    Spts : (ns, 7) int64 cap-point coordinates.
    codim: m - dim(V); number of buckets is 3**codim.
    j    : how many parallel cosets the family unions (e.g. 2 for weight-18, 1 for a
           single 3-flat, 3 for a full 3-flat as 3 parallel 2-flats).

    Returns out (Nsub,) int64: out[s] = sum of the j largest coset occupancies of S
    under direction s.  rigorous max|supp cap S| for the family = out.max().
    """
    nsub = Hs.shape[0]
    ns = Spts.shape[0]
    nb = 1
    for _ in range(codim):
        nb *= P
    out = np.zeros(nsub, dtype=np.int64)
    for s in prange(nsub):
        hist = np.zeros(nb, dtype=np.int64)
        H = Hs[s]
        for p in range(ns):
            cid = 0
            mul = 1
            for k in range(codim):
                acc = 0
                for c in range(7):
                    acc += H[k, c] * Spts[p, c]
                cid += (acc % P) * mul
                mul *= P
            hist[cid] += 1
        # top-j over the nb buckets (j small)
        top = np.zeros(j, dtype=np.int64)   # ascending, top[0] is smallest of the kept
        for b in range(nb):
            v = hist[b]
            if v > top[0]:
                top[0] = v
                # bubble up to keep ascending
                for t in range(1, j):
                    if top[t - 1] > top[t]:
                        tmp = top[t - 1]
                        top[t - 1] = top[t]
                        top[t] = tmp
                    else:
                        break
        ssum = 0
        for t in range(j):
            ssum += top[t]
        out[s] = ssum
    return out


@njit(cache=True, parallel=True)
def inter_counts(supports, mask):
    """|supp cap S| for explicit supports.  supports:(Nsupp,wsupp) point ids; mask:
    (2187,) uint8 indicator of S.  Returns (Nsupp,) int64.  Cross-check only."""
    nsupp, wsupp = supports.shape
    out = np.zeros(nsupp, dtype=np.int64)
    for i in prange(nsupp):
        c = 0
        for t in range(wsupp):
            c += mask[supports[i, t]]
        out[i] = c
    return out
