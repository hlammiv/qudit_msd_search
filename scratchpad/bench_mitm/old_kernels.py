"""Frozen copies of the ORIGINAL (commit 5217ef8) right-hand scan kernels of qmsd.mindist_nb,
with a leading-column range [i_lo, i_hi) added to the ones that lacked it, so the "before"
throughput can be timed on a bounded slice of a production-shaped instance.  The loop bodies are
byte-for-byte the originals (per-probe np.searchsorted into the sorted left table)."""
import numpy as np
from numba import njit, prange

from qmsd.mindist_nb import _fp_dependent, _fp_dependent_b3


@njit(cache=True, parallel=True)
def scan_b3_old(H, p, powers, lcodes, lsupp, i_lo, i_hi):
    """= _scan_right_b3_parallel restricted to leading columns [i_lo, i_hi); returns #verified."""
    r, n = H.shape
    N = lcodes.shape[0]
    a = lsupp.shape[1]
    total = 0
    for i in prange(i_lo, i_hi):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                for c2 in range(1, p):
                    for c3 in range(1, p):
                        target = np.uint64(0)
                        for t in range(r):
                            v = (H[t, i] + c2 * H[t, j] + c3 * H[t, k]) % p
                            target = target + powers[t] * np.uint64((p - v) % p)
                        ii = np.searchsorted(lcodes, target)
                        while ii < N and lcodes[ii] == target:
                            ok = True
                            for s in range(a):
                                ls = lsupp[ii, s]
                                if ls == i or ls == j or ls == k:
                                    ok = False
                                    break
                            if ok and _fp_dependent_b3(H, p, lsupp[ii], i, j, k, a):
                                total += 1
                            ii += 1
    return total


@njit(cache=True, parallel=True)
def scan3_full_old(H, p, powers, lcodes, lsupp, found, witness, i_lo, i_hi):
    """= _scan_right3_full restricted to leading columns [i_lo, i_hi)."""
    r, n = H.shape
    N = lcodes.shape[0]
    for i in prange(i_lo, i_hi):
        rcols = np.empty(3, dtype=np.int32)
        for j in range(i + 1, n):
            if found[0] != 0:
                break
            for k in range(j + 1, n):
                for c1 in range(1, p):
                    for c2 in range(1, p):
                        for c3 in range(1, p):
                            target = np.uint64(0)
                            for t in range(r):
                                v = (c1 * H[t, i] + c2 * H[t, j] + c3 * H[t, k]) % p
                                target = target + powers[t] * np.uint64((p - v) % p)
                            ii = np.searchsorted(lcodes, target)
                            while ii < N and lcodes[ii] == target:
                                l1 = lsupp[ii, 0]; l2 = lsupp[ii, 1]; l3 = lsupp[ii, 2]
                                if (l1 != i and l1 != j and l1 != k and l2 != i and l2 != j
                                        and l2 != k and l3 != i and l3 != j and l3 != k):
                                    rcols[0] = i; rcols[1] = j; rcols[2] = k
                                    if _fp_dependent(H, p, lsupp[ii], rcols, 3, 3):
                                        found[0] = 1
                                        witness[0] = l1; witness[1] = l2; witness[2] = l3
                                        witness[3] = i; witness[4] = j; witness[5] = k
                                ii += 1
    return found[0] != 0


@njit(cache=True, parallel=True)
def census6_old(H, H1, p, powers, lcodes, lsupp, lcoef, i_lo, i_hi):
    """= _census_weight6_3x3 restricted to leading columns [i_lo, i_hi)."""
    r, n = H.shape
    k = H1.shape[0]
    N = lcodes.shape[0]
    fails = np.zeros((n, k), dtype=np.int64)
    found = np.zeros(n, dtype=np.int64)
    for i in prange(i_lo, i_hi):
        for j in range(i + 1, n):
            for kk in range(j + 1, n):
                for c1 in range(1, p):
                    for c2 in range(1, p):
                        for c3 in range(1, p):
                            target = np.uint64(0)
                            for t in range(r):
                                v = (c1 * H[t, i] + c2 * H[t, j] + c3 * H[t, kk]) % p
                                target = target + powers[t] * np.uint64((p - v) % p)
                            ii = np.searchsorted(lcodes, target)
                            while ii < N and lcodes[ii] == target:
                                l1 = lsupp[ii, 0]; l2 = lsupp[ii, 1]; l3 = lsupp[ii, 2]
                                if (l1 != i and l1 != j and l1 != kk and l2 != i and l2 != j
                                        and l2 != kk and l3 != i and l3 != j and l3 != kk):
                                    a2 = lcoef[ii, 0]; a3 = lcoef[ii, 1]
                                    ok = True
                                    for t in range(r):
                                        if (H[t, l1] + a2 * H[t, l2] + a3 * H[t, l3] + c1 * H[t, i]
                                                + c2 * H[t, j] + c3 * H[t, kk]) % p != 0:
                                            ok = False
                                            break
                                    if ok:
                                        found[i] += 1
                                        for o in range(k):
                                            if (H1[o, l1] + a2 * H1[o, l2] + a3 * H1[o, l3]
                                                    + c1 * H1[o, i] + c2 * H1[o, j]
                                                    + c3 * H1[o, kk]) % p != 0:
                                                fails[i, o] += 1
                                ii += 1
    return fails, found


@njit(cache=True, parallel=True)
def census5_old(H, H1, p, powers, lcodes, lsupp, lcoef, i_lo, i_hi):
    """= _census_weight5_per_output restricted to leading columns [i_lo, i_hi)."""
    r, n = H.shape
    k = H1.shape[0]
    N = lcodes.shape[0]
    fails = np.zeros((n, k), dtype=np.int64)
    found = np.zeros(n, dtype=np.int64)
    for i in prange(i_lo, i_hi):
        for j in range(i + 1, n):
            for kk in range(j + 1, n):
                for c2 in range(1, p):
                    for c3 in range(1, p):
                        target = np.uint64(0)
                        for t in range(r):
                            v = (H[t, i] + c2 * H[t, j] + c3 * H[t, kk]) % p
                            target = target + powers[t] * np.uint64((p - v) % p)
                        ii = np.searchsorted(lcodes, target)
                        while ii < N and lcodes[ii] == target:
                            l1 = lsupp[ii, 0]
                            l2 = lsupp[ii, 1]
                            if (l1 != i and l1 != j and l1 != kk and l2 != i and l2 != j
                                    and l2 != kk):
                                a1 = lcoef[ii, 0]
                                a2 = lcoef[ii, 1]
                                ok = True
                                for t in range(r):
                                    if (a1 * H[t, l1] + a2 * H[t, l2] + H[t, i] + c2 * H[t, j]
                                            + c3 * H[t, kk]) % p != 0:
                                        ok = False
                                        break
                                if ok:
                                    found[i] += 1
                                    for o in range(k):
                                        if (a1 * H1[o, l1] + a2 * H1[o, l2] + H1[o, i]
                                                + c2 * H1[o, j] + c3 * H1[o, kk]) % p != 0:
                                            fails[i, o] += 1
                            ii += 1
    return fails, found
