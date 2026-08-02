"""Numba-accelerated meet-in-the-middle weight-d existence kernel for qmsd.mindist.

Same algorithm as the numpy ``_weight_d_exists`` (split d = a + b, halves of size <= 3, match
left syndromes against negated right syndromes by a 64-bit polynomial hash, verify every hash
collision on disjoint supports by an exact F_p rank test), but fused into nopython loops:

  * no ``(r, N)`` syndrome intermediates -- each column-subset's hash is accumulated in a scalar
    loop over the r rows, killing the allocation churn that dominated ``_build_left``/``_encode``;
  * left-table supports stored as int32, not int64 -- half the memory of the giant left table
    (the p=7 wall: (p-1)^a * C(n,a) entries), so the large blocks fit a smaller RAM budget.

The kernel is SELF-CONTAINED: its hash only has to be consistent between the left build and the
right scan (both use ``_POW``), because a collision is never trusted -- it is confirmed by
``_fp_dependent`` before counting. Correctness therefore reduces to "matches brute force", which
is fuzzed in tests. Falls back to the numpy path (imported lazily by the caller) if numba is absent.
"""
from __future__ import annotations

import numpy as np

try:
    from numba import njit
    HAVE_NUMBA = True
except Exception:                       # pragma: no cover - environment without numba
    HAVE_NUMBA = False

    def njit(*args, **kwargs):          # no-op fallback: functions stay pure-python (correct, slow)
        if args and callable(args[0]):
            return args[0]

        def deco(f):
            return f
        return deco


_HASH_BASE = 0x9E3779B97F4A7C15         # same golden-ratio multiplier as qmsd.mindist


def _powers_u64(p, r):
    """Polynomial-hash basis [B^0 .. B^(r-1)] mod 2**64 as uint64 (p unused; signature parity)."""
    out = np.empty(r, dtype=np.uint64)
    mask = (1 << 64) - 1
    v = 1
    for t in range(r):
        out[t] = np.uint64(v)
        v = (v * _HASH_BASE) & mask
    return out


@njit(cache=True)
def _inv(x, p):
    """Modular inverse of x mod p (p prime, small) by trial -- rank-verify inner use only."""
    x %= p
    for y in range(1, p):
        if (x * y) % p == 1:
            return y
    return 0


@njit(cache=True)
def _fp_dependent(H, p, lcols, rcols, a, b):
    """True iff the a+b columns (lcols[:a] then rcols[:b]) of H are F_p-linearly dependent
    (rank < a+b). Exact -- rejects hash-collision false positives, mirrors mindist._dependent."""
    r = H.shape[0]
    d = a + b
    M = np.empty((r, d), dtype=np.int64)
    for s in range(a):
        c = lcols[s]
        for t in range(r):
            M[t, s] = H[t, c] % p
    for s in range(b):
        c = rcols[s]
        for t in range(r):
            M[t, a + s] = H[t, c] % p
    rank = 0
    for col in range(d):
        piv = -1
        for row in range(rank, r):
            if M[row, col] % p != 0:
                piv = row
                break
        if piv < 0:
            continue
        if piv != rank:
            for cc in range(d):
                tmp = M[rank, cc]
                M[rank, cc] = M[piv, cc]
                M[piv, cc] = tmp
        inv = _inv(M[rank, col], p)
        for cc in range(d):
            M[rank, cc] = (M[rank, cc] * inv) % p
        for row in range(r):
            if row != rank and M[row, col] % p != 0:
                f = M[row, col]
                for cc in range(d):
                    val = (M[row, cc] - f * M[rank, cc]) % p
                    if val < 0:
                        val += p
                    M[row, cc] = val
        rank += 1
        if rank == r:
            break
    return rank < d


@njit(cache=True)
def _build_left(H, p, a, powers):
    """All a-subsets (sorted cols) x all nonzero coeff tuples -> (codes_sorted, supp_sorted).
    supp is int32[N, a] (compact); codes is uint64[N]. Fused hash, no (r,N) intermediate."""
    r, n = H.shape
    if a == 1:
        nsub = n
    elif a == 2:
        nsub = n * (n - 1) // 2
    else:
        nsub = n * (n - 1) * (n - 2) // 6
    ncoef = (p - 1) ** a
    N = nsub * ncoef
    codes = np.empty(N, dtype=np.uint64)
    supp = np.empty((N, a), dtype=np.int32)
    idx = 0
    if a == 1:
        for j in range(n):
            for c1 in range(1, p):
                h = np.uint64(0)
                for t in range(r):
                    v = (c1 * H[t, j]) % p
                    h = h + powers[t] * np.uint64(v)
                codes[idx] = h
                supp[idx, 0] = j
                idx += 1
    elif a == 2:
        for j in range(n):
            for k in range(j + 1, n):
                for c1 in range(1, p):
                    for c2 in range(1, p):
                        h = np.uint64(0)
                        for t in range(r):
                            v = (c1 * H[t, j] + c2 * H[t, k]) % p
                            h = h + powers[t] * np.uint64(v)
                        codes[idx] = h
                        supp[idx, 0] = j
                        supp[idx, 1] = k
                        idx += 1
    else:
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    for c1 in range(1, p):
                        for c2 in range(1, p):
                            for c3 in range(1, p):
                                h = np.uint64(0)
                                for t in range(r):
                                    v = (c1 * H[t, i] + c2 * H[t, j] + c3 * H[t, k]) % p
                                    h = h + powers[t] * np.uint64(v)
                                codes[idx] = h
                                supp[idx, 0] = i
                                supp[idx, 1] = j
                                supp[idx, 2] = k
                                idx += 1
    order = np.argsort(codes)
    return codes[order], supp[order]


@njit(cache=True)
def _probe(H, p, lcodes, lsupp, target, rcols, b):
    """Any left entry with code==target whose support is disjoint from rcols[:b] and forms a
    genuine (rank-verified) dependency? Assumes lcodes ascending (searchsorted equal-range)."""
    N = lcodes.shape[0]
    a = lsupp.shape[1]
    ii = np.searchsorted(lcodes, target)
    while ii < N and lcodes[ii] == target:
        disjoint = True
        for s in range(a):
            ls = lsupp[ii, s]
            for u in range(b):
                if ls == rcols[u]:
                    disjoint = False
                    break
            if not disjoint:
                break
        if disjoint and _fp_dependent(H, p, lsupp[ii], rcols, a, b):
            return True
        ii += 1
    return False


@njit(cache=True)
def _scan_right(H, p, b, powers, lcodes, lsupp):
    """True iff some b-subset (first coeff fixed to 1) negated-syndrome hits the left table on a
    disjoint, rank-verified support -- i.e. a genuine weight-(a+b) F_p-dependent column set."""
    r, n = H.shape
    rcols = np.empty(3, dtype=np.int32)
    if b == 1:
        for j in range(n):
            target = np.uint64(0)
            for t in range(r):
                v = H[t, j] % p
                target = target + powers[t] * np.uint64((p - v) % p)
            rcols[0] = j
            if _probe(H, p, lcodes, lsupp, target, rcols, 1):
                return True
    elif b == 2:
        for j in range(n):
            for k in range(j + 1, n):
                for c2 in range(1, p):
                    target = np.uint64(0)
                    for t in range(r):
                        v = (H[t, j] + c2 * H[t, k]) % p
                        target = target + powers[t] * np.uint64((p - v) % p)
                    rcols[0] = j
                    rcols[1] = k
                    if _probe(H, p, lcodes, lsupp, target, rcols, 2):
                        return True
    else:
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    for c2 in range(1, p):
                        for c3 in range(1, p):
                            target = np.uint64(0)
                            for t in range(r):
                                v = (H[t, i] + c2 * H[t, j] + c3 * H[t, k]) % p
                                target = target + powers[t] * np.uint64((p - v) % p)
                            rcols[0] = i
                            rcols[1] = j
                            rcols[2] = k
                            if _probe(H, p, lcodes, lsupp, target, rcols, 3):
                                return True
    return False


def weight_d_exists_fast(H, p, d):
    """True iff H has d columns with a nontrivial zero F_p-combination using all d (assumes no
    weight < d exists -- call with increasing d). Numba MITM; identical result to the numpy path."""
    H = np.ascontiguousarray(np.asarray(H, dtype=np.int64) % p)
    r, n = H.shape
    if d > n:
        return False
    a, b = d // 2, d - d // 2
    powers = _powers_u64(p, r)
    lcodes, lsupp = _build_left(H, p, a, powers)
    if lcodes.shape[0] == 0:
        return False
    return bool(_scan_right(H, p, b, powers, lcodes, lsupp))
