"""Enumerate the affine-flat *directions* (linear subspaces of F_3^7) used by the
structured codeword families, as canonical RREF bases and as syndrome (coset) maps.

A structured low-weight codeword of RM_3(9,7) has support that is a union of cosets
of a common linear *direction* subspace V (dim w).  |supp(c) cap S| is additive over
those cosets, so the rigorous extremum over a whole family reduces to, per direction:
"partition S into the cosets of V and take the top-j coset-occupancy sum".

For that we need, per direction V (dim w):
  - a canonical enumeration of every such V (no double counting), and
  - a syndrome matrix H (codim x 7), codim = m-w, whose row-combination ``H x`` is
    constant exactly on the cosets of V (i.e. ker H = V, rank H = codim).  Then
    coset_id(x) = sum_k (H x)_k * 3^k  in [0, 3^codim).

Counts (Gaussian binomials over F_3^7):  dim 2 -> 99463,  dim 3 -> 925771.
"""
from __future__ import annotations

from itertools import combinations, product

import numpy as np

P = 3
M = 7


def gaussian_binomial(n: int, k: int, q: int = P) -> int:
    num = den = 1
    for i in range(k):
        num *= q ** (n - i) - 1
        den *= q ** (k - i) - 1
    return num // den


def _rref_bases(w: int):
    """Yield every w-dim subspace of F_3^7 exactly once, as a w x 7 RREF basis (int8).

    Canonical form: pick pivot columns p_0<...<p_{w-1}; row i has a 1 at p_i, 0 at the
    other pivot columns and at all columns < p_i, and free F_3 entries at the non-pivot
    columns > p_i.  Each distinct subspace has a unique such matrix.
    """
    cols = range(M)
    for pivots in combinations(cols, w):
        pivset = set(pivots)
        # free (row, col) slots, in a fixed order
        free_slots = []
        for i, pc in enumerate(pivots):
            for c in range(pc + 1, M):
                if c not in pivset:
                    free_slots.append((i, c))
        base = np.zeros((w, M), dtype=np.int8)
        for i, pc in enumerate(pivots):
            base[i, pc] = 1
        if not free_slots:
            yield base.copy()
            continue
        for vals in product(range(P), repeat=len(free_slots)):
            B = base.copy()
            for (i, c), v in zip(free_slots, vals):
                B[i, c] = v
            yield B


def _nullspace_mod3(B: np.ndarray) -> np.ndarray:
    """Basis (codim x 7) of {y : B y = 0 mod 3} for a w x 7 RREF matrix B.

    B is already in RREF, so its pivot columns are the first nonzero of each row; the
    free columns give the standard nullspace basis.  Returns int64, rank = 7 - rank(B).
    """
    B = B.astype(np.int64) % P
    w, n = B.shape
    piv_col = [int(np.nonzero(B[i])[0][0]) for i in range(w)]
    pivset = set(piv_col)
    free = [c for c in range(n) if c not in pivset]
    H = np.zeros((len(free), n), dtype=np.int64)
    for r, f in enumerate(free):
        H[r, f] = 1
        for i, pc in enumerate(piv_col):
            H[r, pc] = (-B[i, f]) % P
    return H


def syndrome_maps(w: int) -> np.ndarray:
    """All direction subspaces of dimension w, as syndrome matrices.

    Returns array (Nsub, codim, 7) int64 with codim = m - w; row block i is H_i with
    ker H_i = V_i.  Nsub = gaussian_binomial(7, w).
    """
    codim = M - w
    Hs = np.empty((gaussian_binomial(M, w), codim, M), dtype=np.int64)
    for i, B in enumerate(_rref_bases(w)):
        Hs[i] = _nullspace_mod3(B)
    return Hs
