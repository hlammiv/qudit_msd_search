"""Memory-rebalanced weight-``d`` codeword counter (additive sibling of qmsd.weightcount).

WHY THIS MODULE EXISTS
----------------------
``qmsd.weightcount.count_weight_d`` counts the weight-``d`` codewords of ``ker(H)`` over F_p
(``= B_d`` of ``G0^perp`` when ``H = G0``) by a meet-in-the-middle split ``d = d1 + d2`` with
``d1 = d//2``.  It builds and SORTS the full weight-``d1`` LEFT table whose every column gets
the FULL nonzero coefficient freedom ``(p-1)**d1``:

    C(n, d1) * (p-1)**d1   entries.

For ``d = 6`` (``d1 = 3``) at the headline pin ``[[237,52,6]]_17`` (p=17, n=237) that is

    C(237,3) * 16**3 = 2_190_670 * 4096 = 8.97e9 entries  ->  ~72 GB of hashes alone,
    far more with supports/coeffs  ->  OOM.

(It is fine for ``d <= 5`` where ``d1 <= 2`` and the LEFT table is small.)

THE REBALANCE (the count analogue of ``qmsd.mindist_balanced``)
--------------------------------------------------------------
We FIX the LEFT block's LEADING (smallest-index) column coefficient to 1 -- i.e. the LEFT
``d1``-subset ranges over coefficients ``(1, c_2, ..., c_{d1})`` with ``c_t`` free nonzero --
so the LEFT table shrinks by a factor ``(p-1)``:

    C(n, d1) * (p-1)**(d1-1) = C(237,3) * 16**2 = 5.61e8 entries  (~4.5 GB hashes; fits),

while the RIGHT half keeps the FULL ``(p-1)**d2`` freedom and is STREAMED first-column by
first-column exactly as ``count_weight_d`` already streams its right (so the ~9e9 right is
never materialised).

THE MULTIPLICITY (derived, then VALIDATED -- see tests/test_weightcount_balanced.py)
------------------------------------------------------------------------------------
``count_weight_d`` finds each weight-``d`` codeword exactly ``C(d,d1)`` times (once per
``d1``/``d2`` support split), so it returns ``raw / C(d,d1) = B_d``.

Fixing the LEFT leading coefficient to 1 changes the multiplicity.  A raw hit is now a pair
``(cw, S_L)`` where ``cw`` is a weight-``d`` codeword, ``S_L`` is a ``d1``-subset of its
support, ``S_R = supp(cw) \\ S_L`` is the (disjoint) ``d2``-subset, AND ``cw`` has coefficient
exactly 1 on ``min(S_L)``.  Group the codewords by their F_p* scaling orbit
``O = {lambda*cw0 : lambda in F_p*}`` (size ``p-1``; every nonzero coordinate of a weight-``d``
word is nonzero, so for each ``d1``-subset ``S_L`` there is exactly ONE ``lambda`` with
``(lambda*cw0)_{min(S_L)} = 1``).  Summing the raw hits over a whole orbit therefore gives

    sum_{lambda} #{S_L : (lambda*cw0)_{min S_L}=1}
      = #{(lambda, S_L) : lambda*(cw0)_{min S_L}=1} = #{d1-subsets S_L} = C(d, d1),

INDEPENDENT of the orbit.  Hence ``raw = C(d,d1) * N_orbits`` and, since
``B_d = (p-1) * N_orbits`` (each orbit is ``p-1`` distinct codewords),

    B_d = (p-1) * raw / C(d, d1).

(Equivalently ``raw_balanced = raw_full / (p-1)``: one fewer coefficient degree of freedom.)
This is the EXACT same ``B_d`` that ``count_weight_d`` and the MacWilliams transform return;
the tests confirm equality on small codes (several p, n, d) where those are feasible.

Like ``count_weight_d`` this is overflow-proof (64-bit polynomial hash of the syndrome, with
the true parity ``H @ cw == 0`` re-verified on every hash hit) and ``dim(G0)``-independent.

This module is ADDITIVE: it reuses the read-only hash primitives of ``qmsd.weightcount`` but
never modifies them.
"""
from __future__ import annotations

import itertools
from math import comb

import numpy as np

# Reuse the (read-only) overflow-proof hash primitives -- never edit weightcount.
from .weightcount import _base_powers, _hash

__all__ = ["count_weight_d_balanced", "left_table_entries"]

_SUPP_DTYPE = np.int32   # column indices; int32 covers n up to ~2.1e9
_COEF_DTYPE = np.int16   # coefficients live in [0, p) with p small


def left_table_entries(n: int, d1: int, p: int) -> int:
    """Number of rebalanced LEFT entries: ``C(n, d1) * (p-1)**(d1-1)``."""
    return comb(n, d1) * (p - 1) ** (d1 - 1)


def _build_left_fixed(H, d1, p, bp):
    """Sorted rebalanced LEFT table of weight-``d1`` partial syndromes.

    Each ``d1``-subset (sorted column indices) is taken with its LEADING (smallest-index)
    column coefficient fixed to 1 and the remaining ``d1-1`` coefficients free nonzero.
    Returns ``(Lh, Lsup, Lcoef)`` sorted by ``Lh``:

      * ``Lh``    : (N,) uint64 polynomial hash of the partial syndrome (overflow-proof).
      * ``Lsup``  : (N, d1) int32 column indices.
      * ``Lcoef`` : (N, d1) int16 coefficients (the leading column's is always 1).

    ``N = C(n, d1) * (p-1)**(d1-1)``.  Memory is preallocated (no list+concatenate doubling).
    """
    rows, n = H.shape
    nz = list(range(1, p))
    combos = np.array(list(itertools.combinations(range(n), d1)), dtype=np.int64)  # (M, d1)
    M = combos.shape[0]
    N = M * (p - 1) ** (d1 - 1)
    Hsel = H[:, combos]  # (rows, M, d1)

    Lh = np.empty(N, dtype=np.uint64)
    Lsup = np.empty((N, d1), dtype=_SUPP_DTYPE)
    Lcoef = np.empty((N, d1), dtype=_COEF_DTYPE)
    combos32 = combos.astype(_SUPP_DTYPE)

    off = 0
    # leading coefficient fixed to 1; the remaining d1-1 coefficients range over nz.
    for tail_cf in itertools.product(nz, repeat=d1 - 1):
        cfa = np.array((1,) + tail_cf, dtype=np.int64)            # length d1, cfa[0] == 1
        s = (Hsel * cfa[None, None, :]).sum(2) % p                # (rows, M)
        Lh[off:off + M] = _hash(s, bp)
        Lsup[off:off + M] = combos32
        Lcoef[off:off + M] = cfa.astype(_COEF_DTYPE)              # broadcast (d1,) -> (M, d1)
        off += M
    assert off == N, (off, N)

    order = np.argsort(Lh, kind="stable")
    return Lh[order], Lsup[order], Lcoef[order]


def count_weight_d_balanced(H, d: int, p: int, progress: int = 0) -> int:
    """Exact number of weight-``d`` codewords of ``ker(H)`` over F_p, memory-rebalanced.

    Identical RESULT to ``qmsd.weightcount.count_weight_d(H, d, p)`` (and to MacWilliams
    ``B_d``), but the weight-``d1`` LEFT table fixes its leading coefficient to 1 so it is
    ``(p-1)x`` smaller and fits in RAM at the ``[[237,52,6]]_17`` pin (where the full LEFT
    table OOMs).  ``d`` is split ``d1 = d//2`` / ``d2 = d - d1``; supports ``2 <= d <= 6``.
    ``progress`` (>0) prints every ``progress`` first-indices.

    Memory ~ the rebalanced LEFT table ``C(n,d1)*(p-1)**(d1-1)``; the RIGHT half (FULL
    coefficients, ``C(n,d2)*(p-1)**d2``) is streamed.  Overflow-proof and ``dim(H)``-independent.
    """
    assert isinstance(p, int) and p >= 2
    assert d >= 1
    H = np.asarray(H, dtype=np.int64) % p
    rows, n = H.shape
    if d > n or rows == 0:
        # rows == 0: ker H = F_p^n, every weight-d pattern is a codeword.
        return comb(n, d) * (p - 1) ** d if rows == 0 else 0

    bp = _base_powers(rows)
    bp3 = bp[None, :, None].astype(np.uint64)
    d1, d2 = d // 2, d - d // 2
    tail = d2 - 1
    if tail > 2:
        raise NotImplementedError("count_weight_d_balanced supports d <= 6 (d2 <= 3)")
    nz = list(range(1, p))
    c1arr = np.array(nz, dtype=np.int64)

    # --- sorted, rebalanced LEFT table (leading coefficient fixed to 1) ---
    Lh, Lsup, Lcoef = _build_left_fixed(H, d1, p, bp)

    # --- stream the FULL-coefficient RIGHT half (identical to count_weight_d's right) ---
    raw = 0
    for i in range(n - d2 + 1):
        if tail == 2:
            jj, kk = np.triu_indices(n - i - 1, 1)
            J, K = jj + i + 1, kk + i + 1
        elif tail == 1:
            J, K = np.arange(i + 1, n), None
        else:  # tail == 0
            J = K = None
        Hi = H[:, i]
        if tail >= 1:
            if J.size == 0:
                continue
            HJ = H[:, J]
            HK = H[:, K] if tail == 2 else None
        for cf in itertools.product(nz, repeat=tail):
            if tail == 2:
                rpart = (cf[0] * HJ + cf[1] * HK) % p
            elif tail == 1:
                rpart = (cf[0] * HJ) % p
            else:
                rpart = np.zeros((rows, 1), dtype=np.int64)
            sc = (c1arr[:, None, None] * Hi[None, :, None] + rpart[None, :, :]) % p  # (nc1,rows,NR)
            with np.errstate(over="ignore"):
                tgt = (bp3 * ((p - sc) % p).astype(np.uint64)).sum(1)                # (nc1,NR)
            pos = np.searchsorted(Lh, tgt, "left")
            valid = pos < Lh.size
            hit = np.zeros(tgt.shape, bool)
            hit[valid] = Lh[pos[valid]] == tgt[valid]
            ha, hr = np.nonzero(hit)
            for a, ri in zip(ha.tolist(), hr.tolist()):
                ci = nz[a]
                t = tgt[a, ri]
                lo = np.searchsorted(Lh, t, "left")
                hi = np.searchsorted(Lh, t, "right")
                if tail == 2:
                    rs, rc = [i, int(J[ri]), int(K[ri])], [ci, cf[0], cf[1]]
                elif tail == 1:
                    rs, rc = [i, int(J[ri])], [ci, cf[0]]
                else:
                    rs, rc = [i], [ci]
                rset = set(rs)
                for li in range(lo, hi):
                    ls = Lsup[li]
                    if rset.isdisjoint(int(x) for x in ls):
                        cols = np.array([int(x) for x in ls] + rs)
                        coef = np.array(
                            [int(x) for x in Lcoef[li]] + rc, dtype=np.int64
                        )
                        if np.any((H[:, cols] @ coef) % p):
                            continue  # reject hash collision
                        raw += 1
        if progress and i % progress == 0:
            est = (p - 1) * raw / comb(d, d1)
            print(f"  i={i}/{n - d2 + 1} B_d~={est:.0f}", flush=True)

    cdd1 = comb(d, d1)
    assert raw % cdd1 == 0, f"raw={raw} not divisible by C({d},{d1})={cdd1}"
    return (p - 1) * (raw // cdd1)
