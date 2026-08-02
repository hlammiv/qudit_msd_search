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


@njit(cache=True)
def _scan_right_b4(H, p, powers, lcodes, lsupp, i_max):
    """b=4 right stream (first coeff fixed to 1) against an a=2 left table -> weight-6 witness.
    This is the 2+4 split for d=6: only the a=2 left table is stored (fits), the weight-4 side is
    streamed with early witness-exit, so d=6 is certifiable WITHOUT the a=3 table (~TBs at n~2000).
    Full scan is C(n,4)(p-1)^3 (proving d>=7 stays infeasible); finding a weight-6 witness exits
    early. ``i_max`` bounds the leading column so a d>=7 set (no witness) can't run forever -- a
    d=6 code has weight-6 words with small leading index and is found in a modest budget. Assumes
    no weight<6 exists (call after the d<=5 search returns nothing)."""
    r, n = H.shape
    rcols = np.empty(4, dtype=np.int32)
    for i in range(min(i_max, n)):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                for l in range(k + 1, n):
                    for c2 in range(1, p):
                        for c3 in range(1, p):
                            for c4 in range(1, p):
                                target = np.uint64(0)
                                for t in range(r):
                                    v = (H[t, i] + c2 * H[t, j] + c3 * H[t, k] + c4 * H[t, l]) % p
                                    target = target + powers[t] * np.uint64((p - v) % p)
                                rcols[0] = i; rcols[1] = j; rcols[2] = k; rcols[3] = l
                                if _probe(H, p, lcodes, lsupp, target, rcols, 4):
                                    return True
    return False


def weight6_exists_via_2_4(H, p, i_max=None):
    """True iff H has 6 F_p-dependent columns, found via a 2+4 split (a=2 left table + b=4 stream).
    Memory = the a=2 left table only, so this certifies d=6 at n~2000 where the standard a=3,b=3
    d=6 path OOMs. ``i_max`` bounds the weight-4 leading column (default: all n) so a d>=7 set
    can't hang; a genuine d=6 code is found within a small budget. Assumes no weight<6 (call after
    min_dependent_columns(d_max=5) finds nothing)."""
    H = np.ascontiguousarray(np.asarray(H, dtype=np.int64) % p)
    r, n = H.shape
    if n < 6:
        return False
    powers = _powers_u64(p, r)
    lcodes, lsupp = _build_left(H, p, 2, powers)
    if lcodes.shape[0] == 0:
        return False
    return bool(_scan_right_b4(H, p, powers, lcodes, lsupp, n if i_max is None else int(i_max)))


def min_distance_upto6_lowmem(H, p):
    """Exact min distance as an int in 1..6, or 7 meaning ">=7", using ONLY a<=2 left tables.

    d=1..5 go through the standard MITM halves (a<=2, memory-safe); d=6 uses the 2+4 split
    (a=2 left table + b=4 stream) instead of the a=3,b=3 path whose C(n,3)(p-1)^3 left table is
    ~TBs at n~2000. So this certifies d=6 codes on a normal RAM budget where min_dependent_columns
    OOMs -- the enabler for the gamma<1 (d>=6) regime. Returning 7 requires the full b=4 scan
    (feasible only at small n); at large n a d>=7 code makes the b=4 stream infeasible (expected)."""
    H = np.ascontiguousarray(np.asarray(H, dtype=np.int64) % p)
    r, n = H.shape
    if r == 0:
        return 1
    if (H == 0).all(axis=0).any():
        return 1
    for d in (2, 3, 4, 5):
        if weight_d_exists_fast(H, p, d):
            return d
    if weight6_exists_via_2_4(H, p):
        return 6
    return 7


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
