"""Sort + merge-join alternative for the 2+3 right scan (bench-only, for comparison with the
bucket-indexed kernel that went into qmsd.mindist_nb).

Per leading column (prange): the right-hand targets are generated into a buffer of ``cap``
entries (hash + packed (j, k, c2, c3)); when full the buffer is argsorted (numba quicksort) and
joined against the ASCENDING left table with a galloping merge (both sides touched in increasing
address order), each hash hit verified exactly as in the production kernels.  Same probe SET as
``_scan_right_b3_idx`` (the order inside a buffer differs, which is irrelevant for a count).

Memory per thread: cap * (8 + 8 + 8) B (keys, meta, argsort permutation) -> 24 MB at cap = 2^20.
"""
import numpy as np
from numba import njit, prange

from qmsd.mindist_nb import _fp_dependent_b3, _scaled_cols, _mod_tables


@njit(cache=True)
def _flush_merge(H, p, lcodes, lsupp, keys, meta, m, i):
    """Join keys[:m] (unsorted) against lcodes (sorted) and return the verified-hit count."""
    N = lcodes.shape[0]
    a = lsupp.shape[1]
    order = np.argsort(keys[:m])
    hits = 0
    pos = 0
    for q in range(m):
        o = order[q]
        target = keys[o]
        # gallop from pos: find first idx >= pos with lcodes[idx] >= target
        if pos < N and lcodes[pos] < target:
            step = 1
            lo = pos
            while pos + step < N and lcodes[pos + step] < target:
                lo = pos + step
                step *= 2
            hi = pos + step
            if hi > N:
                hi = N
            # binary search in (lo, hi]: lcodes[lo] < target <= lcodes[hi] (or hi == N)
            while hi - lo > 1:
                mid = (lo + hi) >> 1
                if lcodes[mid] < target:
                    lo = mid
                else:
                    hi = mid
            pos = hi
        ii = pos
        if ii < N and lcodes[ii] == target:
            mm = meta[o]
            j = np.int64(mm >> np.int64(40))
            k = np.int64((mm >> np.int64(16)) & np.int64(0xFFFFFF))
            while ii < N and lcodes[ii] == target:
                ok = True
                for s in range(a):
                    ls = lsupp[ii, s]
                    if ls == i or ls == j or ls == k:
                        ok = False
                        break
                if ok and _fp_dependent_b3(H, p, lsupp[ii], i, j, k, a):
                    hits += 1
                ii += 1
    return hits


@njit(cache=True, parallel=True)
def scan_b3_merge(H, p, powers, lcodes, lsupp, i_lo, i_hi, cap):
    """Sort+merge 2+3 scan over leading columns [i_lo, i_hi): number of verified hits."""
    r, n = H.shape
    CHT = _scaled_cols(H, p)
    red, neg = _mod_tables(p)
    total = 0
    for i in prange(i_lo, min(i_hi, n)):
        U = np.empty((p, r), dtype=np.int64)
        keys = np.empty(cap, dtype=np.uint64)
        meta = np.empty(cap, dtype=np.int64)
        m = 0
        for j in range(i + 1, n):
            for c2 in range(1, p):
                for t in range(r):
                    U[c2, t] = red[CHT[i, 1, t] + CHT[j, c2, t]]
            for k in range(j + 1, n):
                for c2 in range(1, p):
                    for c3 in range(1, p):
                        target = np.uint64(0)
                        for t in range(r):
                            target = target + powers[t] * np.uint64(neg[U[c2, t] + CHT[k, c3, t]])
                        keys[m] = target
                        meta[m] = (np.int64(j) << np.int64(40)) | (np.int64(k) << np.int64(16)) \
                            | (np.int64(c2) << np.int64(8)) | np.int64(c3)
                        m += 1
                        if m == cap:
                            total += _flush_merge(H, p, lcodes, lsupp, keys, meta, m, i)
                            m = 0
        if m > 0:
            total += _flush_merge(H, p, lcodes, lsupp, keys, meta, m, i)
    return total
