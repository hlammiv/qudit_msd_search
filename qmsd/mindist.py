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
vector in F_p^r) is encoded as a single int64 by mixed radix, which is a *bijection*, so a
code match is an exact syndrome match -- no false positives, only a disjoint-support test.

This is exact and certified: it returns the true minimum distance (up to d_max), never an
under-report. Supports d up to 6 (halves of size <= 3), which covers every code in the
paper's Table 3. Pure numpy; no external dependencies.
"""
from __future__ import annotations

import itertools

import numpy as np

from .triorthogonal import dual_matrix

# int64 must hold p**r for the syndrome encoding to be injective.
_INT64_MAX = np.iinfo(np.int64).max


def _as_int_H(generator, p):
    """Parity check H = dual(generator) as an int ndarray, reduced mod p."""
    H = np.asarray(dual_matrix(generator, p)) % p
    return H.astype(np.int64, copy=False)


def _powers(p, r):
    if r == 0:
        return np.ones(0, dtype=np.int64)
    if p ** r > _INT64_MAX:
        raise OverflowError(
            f"syndrome encoding p**r = {p}**{r} exceeds int64; redundancy too large "
            f"for the MITM encoder"
        )
    return (p ** np.arange(r, dtype=np.int64)).astype(np.int64)


def _encode(S, powers):
    """Encode syndrome columns S (shape (r, N), entries in [0,p)) to int64 codes (N,)."""
    if S.shape[0] == 0:
        return np.zeros(S.shape[1], dtype=np.int64)
    return powers @ S  # (r,)·(r,N) -> (N,)


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

    codes = np.concatenate(code_chunks) if code_chunks else np.empty(0, np.int64)
    supps = np.concatenate(supp_chunks) if supp_chunks else np.empty((0, a), np.int64)
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
                    return True
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
    powers = _powers(p, r)
    # Search weight 1..cap. Any r+1 columns of an r-row matrix are dependent, so the true
    # distance is at most r+1; the MITM halves (each size <= 3) only certify up to HARD_CAP,
    # so cap the search there and raise (never mis-report) if nothing is found by then.
    HARD_CAP = 6
    cap = min(d_max if d_max is not None else (r + 1), n, HARD_CAP)
    for d in range(1, cap + 1):
        if d == 1:
            if (H == 0).all(axis=0).any():
                return 1
            continue
        if _weight_d_exists(H, p, d, powers):
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
