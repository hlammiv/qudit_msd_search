"""Meet-in-the-middle exact minimum distance for the punctured-RM dual codes.

The quantum distance of a triorthogonal code is d(G0^perp) (NOTES Thm 3), i.e. the
minimum weight of a nonzero codeword of the code whose PARITY CHECK is G0. Equivalently
it is the minimum number of F_p-linearly-dependent columns of the (small) r x n matrix
H = G0, where r = dim(G0) is the code's redundancy (small for these codes, ~16-38).

A weight-d codeword is a choice of d columns i_1<...<i_d and nonzero coefficients
a_1,...,a_d in F_p with  sum_t a_t H[:,i_t] = 0.  Naively testing this scans ~C(n,d)
column subsets. Meet-in-the-middle splits d = a + b:

    s_L = sum over an a-subset (free nonzero coeffs)            -> hash
    s_R = sum over a disjoint b-subset (first coeff fixed to 1) -> stream

and a codeword exists iff some s_L = -s_R on disjoint supports. Fixing the b-block's
leading coefficient to 1 removes the global scaling redundancy while still covering every
codeword (scale it so its b-block's smallest column has coefficient 1). Each syndrome (a
vector in F_p^r) is matched by a 64-bit polynomial HASH; a hash collision on disjoint
supports is VERIFIED exactly by an F_p rank test (_dependent) before it counts, so there are
no false positives and the redundancy r is UNBOUNDED (no p-adic int64 overflow -- the wall
we kept hitting at large p / mid-k).

This is exact and certified: it returns the true minimum distance (up to d_max), never an
under-report. Supports d up to 6 (halves of size <= 3), which covers every code in the
paper's Table 3. Pure numpy; no external dependencies.
"""
from __future__ import annotations

import itertools

import numpy as np

from .triorthogonal import dual_matrix
from .mindist_nb import weight_d_exists_fast, HAVE_NUMBA as _HAVE_NUMBA

# Syndromes are matched by a 64-bit polynomial HASH (not the p-adic int64 radix, which
# overflows once p**r > 2**63 -- the wall we kept hitting at large p / large redundancy).
# The hash is NOT injective, so every collision on disjoint supports is VERIFIED exactly by
# an F_p rank test (_dependent) before it counts -- results stay certified, r is unbounded.
_HASH_BASE = np.uint64(0x9E3779B97F4A7C15)  # fixed 64-bit golden-ratio multiplier


def _as_int_H(generator, p):
    """Parity check H = dual(generator) as an int ndarray, reduced mod p."""
    H = np.asarray(dual_matrix(generator, p)) % p
    return H.astype(np.int64, copy=False)


def _powers(p, r):
    """Polynomial-hash basis [B^0 .. B^(r-1)] mod 2**64 (uint64). Overflow-proof replacement
    for the p-adic radix; no cap on the redundancy r (p is unused, kept for signature compat)."""
    if r == 0:
        return np.ones(0, dtype=np.uint64)
    B = int(_HASH_BASE)
    mask = (1 << 64) - 1
    vals = [1]
    for _ in range(1, r):
        vals.append((vals[-1] * B) & mask)
    return np.array(vals, dtype=np.uint64)


def _encode(S, powers):
    """Hash syndrome columns S (shape (r, N), entries in [0,p)) -> (N,) uint64 (wraps mod 2**64)."""
    if S.shape[0] == 0:
        return np.zeros(S.shape[1], dtype=np.uint64)
    with np.errstate(over="ignore"):
        return (powers[:, None] * S.astype(np.uint64)).sum(axis=0)


def _fp_rank(rows_mat, p) -> int:
    """Rank over F_p of a small integer matrix given as rows."""
    M = np.asarray(rows_mat, dtype=np.int64) % p
    nr, nc = M.shape
    rank = 0
    for c in range(nc):
        piv = next((i for i in range(rank, nr) if M[i, c] % p), None)
        if piv is None:
            continue
        M[[rank, piv]] = M[[piv, rank]]
        M[rank] = (M[rank] * pow(int(M[rank, c]), p - 2, p)) % p
        for i in range(nr):
            if i != rank and M[i, c] % p:
                M[i] = (M[i] - M[i, c] * M[rank]) % p
        rank += 1
        if rank == nr:
            break
    return rank


def _dependent(H, p, cols) -> bool:
    """True iff columns H[:, cols] are F_p-linearly dependent (rank < #cols). Given all lower
    weights are already ruled out, this certifies a genuine full-support weight-|cols| codeword
    and rejects hash-collision false positives."""
    cols = list(cols)
    return _fp_rank(H[:, cols].T, p) < len(cols)


def _pairs_after(start, n):
    """All (j,k) with start <= j < k < n, as two index arrays."""
    cols = np.arange(start, n)
    if cols.size < 2:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    j, k = np.triu_indices(cols.size, k=1)
    return cols[j], cols[k]


def _build_left(H, p, a, powers):
    """All a-subsets with FULL free nonzero coeffs -> (codes_sorted, supports_sorted).

    supports is an (N, a) int array of column indices (sorted within a row). The returned
    arrays are sorted by code so a target can be range-looked-up via searchsorted.
    """
    r, n = H.shape
    nz = range(1, p)
    code_chunks, supp_chunks = [], []

    if a == 1:
        cols = np.arange(n)
        for c in nz:
            S = (c * H[:, cols]) % p
            code_chunks.append(_encode(S, powers))
            supp_chunks.append(cols[:, None].copy())
    elif a == 2:
        j, k = _pairs_after(0, n)
        for c1 in nz:
            Hj = (c1 * H[:, j]) % p
            for c2 in nz:
                S = (Hj + c2 * H[:, k]) % p
                code_chunks.append(_encode(S, powers))
                supp_chunks.append(np.stack([j, k], axis=1))
    elif a == 3:
        for i in range(n - 2):
            j, k = _pairs_after(i + 1, n)
            if j.size == 0:
                continue
            Hjk_cache = {}
            for c1 in nz:
                Hi = (c1 * H[:, i]) % p
                for c2 in nz:
                    Hj = (c2 * H[:, j]) % p
                    for c3 in nz:
                        S = (Hi[:, None] + Hj + c3 * H[:, k]) % p
                        code_chunks.append(_encode(S, powers))
                        ii = np.full(j.size, i)
                        supp_chunks.append(np.stack([ii, j, k], axis=1))
    else:
        raise ValueError("left half size must be 1, 2, or 3")

    codes = np.concatenate(code_chunks) if code_chunks else np.empty(0, np.uint64)
    # int32 supports halve the giant left table's memory (columns << 2**31 for any certifiable
    # block); indices are cast to int downstream so the narrower dtype is transparent.
    supps = (np.concatenate(supp_chunks).astype(np.int32, copy=False)
             if supp_chunks else np.empty((0, a), np.int32))
    order = np.argsort(codes, kind="stable")
    return codes[order], supps[order]


def _iter_right(H, p, b, powers):
    """Yield (target_codes, supports) batches for all b-subsets, first coeff fixed to 1.

    ``target_codes`` is encode(-s_R) -- the LEFT code we must find for a collision.
    supports is (batch, b) column indices.
    """
    r, n = H.shape
    nz = range(1, p)

    def neg_code(S):
        return _encode((p - S) % p, powers)

    if b == 1:
        cols = np.arange(n)
        S = H[:, cols] % p
        yield neg_code(S), cols[:, None]
    elif b == 2:
        j, k = _pairs_after(0, n)
        Hj = H[:, j] % p
        supp = np.stack([j, k], axis=1)
        for c2 in nz:
            S = (Hj + c2 * H[:, k]) % p
            yield neg_code(S), supp
    elif b == 3:
        for i in range(n - 2):
            j, k = _pairs_after(i + 1, n)
            if j.size == 0:
                continue
            Hi = H[:, i] % p
            supp = np.stack([np.full(j.size, i), j, k], axis=1)
            for c2 in nz:
                Hj = (c2 * H[:, j]) % p
                for c3 in nz:
                    S = (Hi[:, None] + Hj + c3 * H[:, k]) % p
                    yield neg_code(S), supp
    else:
        raise ValueError("right half size must be 1, 2, or 3")


def _weight_d_exists(H, p, d, powers):
    """True iff H has d columns with a nontrivial zero F_p-combination using all d.

    Assumes no weight < d codeword exists (call with increasing d), so a code/syndrome
    collision on DISJOINT supports is a genuine weight-d codeword.
    """
    r, n = H.shape
    if d > n:
        return False
    a, b = d // 2, d - d // 2  # a <= b ; halves of size <= 3 (so d <= 6)

    left_codes, left_supps = _build_left(H, p, a, powers)
    if left_codes.size == 0:
        return False

    for target_codes, right_supps in _iter_right(H, p, b, powers):
        # Vectorised membership: which targets appear among the sorted left codes?
        pos = np.searchsorted(left_codes, target_codes, side="left")
        in_range = pos < left_codes.size
        hit = np.zeros(target_codes.shape, dtype=bool)
        hit[in_range] = left_codes[pos[in_range]] == target_codes[in_range]
        if not hit.any():
            continue
        # Resolve the (rare) hits: every left entry sharing the code, disjoint-support test.
        for ridx in np.nonzero(hit)[0]:
            t = target_codes[ridx]
            lo = np.searchsorted(left_codes, t, side="left")
            hi = np.searchsorted(left_codes, t, side="right")
            rcols = set(int(x) for x in right_supps[ridx])
            for li in range(lo, hi):
                lcols = left_supps[li]
                if rcols.isdisjoint(int(x) for x in lcols):
                    if _dependent(H, p, [int(x) for x in lcols] + list(rcols)):
                        return True  # exact rank check rejects hash-collision false positives
    return False


def min_dependent_columns(H, p, d_max=None) -> int:
    """Minimum number of F_p-linearly-dependent columns of H (= distance of ker H).

    Exact and certified up to d_max (default: the redundancy bound r+1, since any r+1
    columns of an r-row matrix are dependent). Raises ValueError if the trivial code
    (no dependent set <= d_max) -- which cannot happen for d_max >= r+1.
    """
    H = np.asarray(H, dtype=np.int64) % p
    r, n = H.shape
    if r == 0:
        return 1  # empty parity check: every single column is a codeword
    powers = _powers(p, r) if not _HAVE_NUMBA else None
    # Search weight 1..cap. Any r+1 columns of an r-row matrix are dependent, so the true
    # distance is at most r+1; the MITM halves (each size <= 3) only certify up to HARD_CAP,
    # so cap the search there and raise (never mis-report) if nothing is found by then.
    # The numba kernel (weight_d_exists_fast) fuses the syndrome build+hash into nopython loops
    # (~3x, compact int32 supports) and gives the IDENTICAL result -- the numpy _weight_d_exists
    # is the reference/fallback when numba is unavailable (fuzzed equal in tests).
    HARD_CAP = 6
    cap = min(d_max if d_max is not None else (r + 1), n, HARD_CAP)
    for d in range(1, cap + 1):
        if d == 1:
            if (H == 0).all(axis=0).any():
                return 1
            continue
        found = (weight_d_exists_fast(H, p, d) if _HAVE_NUMBA
                 else _weight_d_exists(H, p, d, powers))
        if found:
            return d
    raise ValueError(
        f"no dependent column set of size <= {cap} found; minimum distance exceeds {cap} "
        f"(MITM certifies up to {HARD_CAP}) -- not certified"
    )


# ---------------------------------------------------------------------------
# Parallel MITM: identical result to min_dependent_columns, but the weight-3 right-half
# enumeration is split across processes by the right subset's first-column index. The
# read-only left table is auto-memmapped by joblib (one shared copy on disk, mmapped by
# every worker) so peak RAM is ~the single-thread footprint, NOT n_jobs copies. Hard-capped
# at weight 3 (d <= 6): it never builds the weight-4 (C(n,4)) table that causes OOM.
# ---------------------------------------------------------------------------
def _right_i_collision(H, p, powers, left_codes, left_supps, i_list):
    """Worker: does any weight-3 right-subset whose first column is in ``i_list`` collide
    with the shared, read-only left table on disjoint supports? Mirrors the b==3 branch of
    _iter_right + the matching in _weight_d_exists, restricted to ``i_list``."""
    r, n = H.shape
    nz = range(1, p)
    for i in i_list:
        j, k = _pairs_after(i + 1, n)
        if j.size == 0:
            continue
        Hi = H[:, i] % p
        supp = np.stack([np.full(j.size, i), j, k], axis=1)
        for c2 in nz:
            Hj = (c2 * H[:, j]) % p
            for c3 in nz:
                S = (Hi[:, None] + Hj + c3 * H[:, k]) % p
                target = _encode((p - S) % p, powers)
                pos = np.searchsorted(left_codes, target, side="left")
                in_range = pos < left_codes.size
                hit = np.zeros(target.shape, dtype=bool)
                hit[in_range] = left_codes[pos[in_range]] == target[in_range]
                if not hit.any():
                    continue
                for ridx in np.nonzero(hit)[0]:
                    t = target[ridx]
                    lo = np.searchsorted(left_codes, t, side="left")
                    hi = np.searchsorted(left_codes, t, side="right")
                    rcols = set(int(x) for x in supp[ridx])
                    for li in range(lo, hi):
                        lcols = left_supps[li]
                        if rcols.isdisjoint(int(x) for x in lcols):
                            if _dependent(H, p, [int(x) for x in lcols] + list(rcols)):
                                return True  # exact rank check rejects hash collisions
    return False


def _weight_d_exists_parallel(H, p, d, powers, n_jobs):
    r, n = H.shape
    if d > n:
        return False
    a, b = d // 2, d - d // 2
    left_codes, left_supps = _build_left(H, p, a, powers)
    if left_codes.size == 0:
        return False
    if b != 3 or n_jobs == 1:
        return _weight_d_exists(H, p, d, powers)  # b<=2 is cheap -> serial
    from joblib import Parallel, delayed
    # round-robin i-slices balance load (small i has the widest (j,k) range -> most work)
    chunks = [list(range(w, n - 2, n_jobs)) for w in range(n_jobs)]
    res = Parallel(n_jobs=n_jobs)(
        delayed(_right_i_collision)(H, p, powers, left_codes, left_supps, ch)
        for ch in chunks if ch
    )
    return any(res)


def min_dependent_columns_parallel(H, p, d_max=None, n_jobs=-1) -> int:
    """Parallel min_dependent_columns: identical result, weight-3 right-half search split
    across ``n_jobs`` processes. joblib auto-memmaps the read-only left table so RAM stays
    ~one shared copy (not n_jobs copies); the weight-3 cap (d<=6) keeps it out of the
    weight-4 memory-blowup regime. n_jobs=-1 uses all cores."""
    import os
    H = np.asarray(H, dtype=np.int64) % p
    r, n = H.shape
    if r == 0:
        return 1
    powers = _powers(p, r)
    HARD_CAP = 6
    cap = min(d_max if d_max is not None else (r + 1), n, HARD_CAP)
    if n_jobs in (-1, None):
        n_jobs = os.cpu_count() or 1
    for d in range(1, cap + 1):
        if d == 1:
            if (H == 0).all(axis=0).any():
                return 1
            continue
        if _weight_d_exists_parallel(H, p, d, powers, n_jobs):
            return d
    raise ValueError(
        f"no dependent column set of size <= {cap} found; minimum distance exceeds {cap} "
        f"(MITM certifies up to {HARD_CAP}) -- not certified"
    )


def min_distance_certified(generator, p, d_max=None) -> int:
    """Exact minimum distance of the code spanned by ``generator``, via MITM on its dual.

    Equivalent to qmsd.distance.min_distance but uses the meet-in-the-middle column-
    dependency search on the (small-redundancy) parity check -- far faster for the
    high-rate punctured-RM duals, so it certifies the large search codes the naive scan
    cannot. Distances up to 6 (which bounds every Table-3 code).
    """
    H = _as_int_H(generator, p)
    return min_dependent_columns(H, p, d_max=d_max)


# ---------------------------------------------------------------------------
# Brute-force reference (ground truth for verification / small codes only).
# ---------------------------------------------------------------------------
def min_distance_bruteforce(generator, p, d_max=None) -> int:
    """Ground-truth minimum distance by exhaustive column-dependency search (slow).

    For every d = 1, 2, ..., test all C(n,d) column subsets for F_p-linear dependence
    (rank < d). O(C(n,d)) -- use only for small codes / verification.
    """
    H = _as_int_H(generator, p)
    r, n = H.shape
    if r == 0:
        return 1
    cap = (r + 1) if d_max is None else min(d_max, n)
    for d in range(1, cap + 1):
        for cols in itertools.combinations(range(n), d):
            sub = H[:, list(cols)]
            if np.linalg.matrix_rank(_GF(p)(sub.astype(int) % p)) < d:
                return d
    raise ValueError(f"distance exceeds d_max={cap}")


def _GF(p):
    import galois
    return galois.GF(p)
