"""Direct, overflow-proof weight-d codeword counter for ``G0^perp`` (exact ``A_d``).

For a punctured triorthogonal code the quantum distance is the minimum weight of the dual
``G0^perp`` and the paper's ``A_d`` counts its minimum-weight (logical) operators.  In the
validated regime -- where every minimum-weight dual codeword is logical (the stabilizer code
``Gp^perp`` has minimum weight ``> d``; empirically certified for the Table-3 codes, see
``weightdist``) -- ``A_d = B_d``, the total number of weight-``d`` codewords of ``G0^perp``.

``exact_distance_and_Ad`` (``qmsd.weightdist``) gets ``B_d`` from the MacWilliams transform of
the *enumerated* generator ``G0``; that costs ``p**dim(G0)`` and is INFEASIBLE once ``dim(G0)``
is large (e.g. the flagship ``[[519,106,5]]_5`` has ``dim G0 = 16`` ⇒ ``5**16 ≈ 1.5e11``).

This module instead counts the weight-``d`` codewords *directly* as the weight-``d``
F_p-dependent column sets of the parity check ``H = G0`` (``ker H = G0^perp``), by a
meet-in-the-middle split ``d = d1 + d2``:

  * build the small SORTED table of all weight-``d1`` partial syndromes (``C(n,d1)(p-1)^{d1}``),
  * STREAM the weight-``d2`` right half by first-column index ``i`` so the (huge) ``C(n,d2)``
    array never materialises at once,
  * match by a 64-bit polynomial HASH of the length-``rows`` syndrome (overflow-proof, unlike the
    ``p**rows`` int64 p-adic), VERIFYING the true parity ``H @ cw == 0`` on every hash hit so a
    collision is simply rejected.

Each weight-``d`` codeword is found exactly ``C(d,d1)`` times (its support's ``d1``/``d2`` splits),
so ``A_d = (#verified matches) / C(d,d1)``.  Cost ``~ C(n,d2)(p-1)^{d2}`` syndrome evaluations,
INDEPENDENT of ``dim(G0)`` -- the regime where MacWilliams is infeasible.
"""
from __future__ import annotations

import itertools
from math import comb

import numpy as np

__all__ = ["count_weight_d"]

_HASH_BASE = np.uint64(0x9E3779B97F4A7C15)  # fixed 64-bit golden-ratio multiplier


def _base_powers(rows: int) -> np.ndarray:
    """``[B^0..B^(rows-1)]`` mod 2**64 (uint64) -- the polynomial-hash basis (explicit wrap)."""
    B = int(_HASH_BASE)
    mask = (1 << 64) - 1
    vals = [1]
    for _ in range(1, rows):
        vals.append((vals[-1] * B) & mask)
    return np.array(vals, dtype=np.uint64)


def _hash(s: np.ndarray, bp: np.ndarray) -> np.ndarray:
    """Hash each column of ``s`` (rows x N, entries in [0,p)) -> (N,) uint64 (wraps mod 2**64)."""
    with np.errstate(over="ignore"):
        return (bp[:, None] * s.astype(np.uint64)).sum(axis=0)


def count_weight_d(H, d: int, p: int, progress: int = 0) -> int:
    """Exact number of weight-``d`` codewords of ``ker(H)`` over F_p (``H`` = ``rows x n``).

    For ``H = G0`` this is ``B_d`` of ``G0^perp`` (= the paper's ``A_d`` in the validated regime).
    Overflow-proof and ``dim(G0)``-independent (no ``p**dim`` enumeration).  ``d`` is split
    ``d1 = d//2`` / ``d2 = d - d1``; supports ``2 <= d``.  ``progress`` (>0) prints every
    ``progress`` first-indices.  Memory ~ the weight-``d1`` table (``C(n,d1)(p-1)^{d1}``); the
    weight-``d2`` half is streamed.
    """
    assert isinstance(p, int) and p >= 2
    assert d >= 1
    H = np.asarray(H, dtype=np.int64) % p
    rows, n = H.shape
    if d > n or rows == 0:
        # rows == 0: ker H = F_p^n, every weight-d pattern is a codeword
        return comb(n, d) * (p - 1) ** d if rows == 0 else 0
    bp = _base_powers(rows)
    bp3 = bp[None, :, None].astype(np.uint64)
    d1, d2 = d // 2, d - d // 2
    tail = d2 - 1
    nz = list(range(1, p))
    c1arr = np.array(nz, dtype=np.int64)

    # --- sorted left table: weight-d1 partial-syndrome hashes ---
    combos = np.array(list(itertools.combinations(range(n), d1)), dtype=np.int64)
    Hsel = H[:, combos]
    Lh, Lsup, Lcoef = [], [], []
    for cf in itertools.product(nz, repeat=d1):
        cfa = np.array(cf, dtype=np.int64)
        s = (Hsel * cfa[None, None, :]).sum(2) % p
        Lh.append(_hash(s, bp))
        Lsup.append(combos)
        Lcoef.append(np.broadcast_to(cfa, combos.shape))
    Lh = np.concatenate(Lh)
    Lsup = np.concatenate(Lsup)
    Lcoef = np.concatenate(Lcoef)
    order = np.argsort(Lh, kind="stable")
    Lh, Lsup, Lcoef = Lh[order], Lsup[order], Lcoef[order]

    raw = 0
    for i in range(n - d2 + 1):
        if tail == 2:
            jj, kk = np.triu_indices(n - i - 1, 1)
            J, K = jj + i + 1, kk + i + 1
        elif tail == 1:
            J, K = np.arange(i + 1, n), None
        elif tail == 0:
            J = K = None
        else:
            raise NotImplementedError("count_weight_d supports d <= 6 (d2 <= 3)")
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
                        cols = np.array(list(ls) + rs)
                        coef = np.array(list(Lcoef[li]) + rc)
                        if np.any((H[:, cols] @ coef) % p):
                            continue  # reject hash collision
                        raw += 1
        if progress and i % progress == 0:
            print(f"  i={i}/{n - d2 + 1} A_d~={raw / comb(d, d1):.0f}", flush=True)

    assert raw % comb(d, d1) == 0, f"raw={raw} not divisible by C({d},{d1})={comb(d, d1)}"
    return raw // comb(d, d1)
