"""Puncture-column <-> point indexing and Manhattan puncture sets (NOTES Appendix C; sec 5).

The unpunctured RM_p(r,m) generator matrix has p^m columns, **1-indexed**
(first column c=1). Appendix C fixes the column<->point bijection: column c
corresponds to the point x = (x_1,...,x_m) in F_p^m given by the base-p
expansion of c-1 with x_1 least significant,

    c - 1 = x_1 + x_2 p + x_3 p^2 + ... + x_m p^(m-1),   x_i in {0,...,p-1},
    x_i   = floor((c-1) / p^(i-1)) mod p.

Puncturing by Manhattan weight uses the set (NOTES sec 5)

    S_p(w) = { x in F_p^m : |x|_M = sum_i x_i <= w },   |S_p(w)| = [m,<=w]_p.
"""
from __future__ import annotations


def column_to_point(c, m, p) -> tuple:
    """1-indexed column c -> point (x_1..x_m) in F_p^m, x_1 least significant (NOTES Appendix C).

    Decodes the base-p digits of c-1: x_i = floor((c-1)/p^(i-1)) mod p.
    """
    assert p >= 2, "p must be >= 2"
    assert m >= 1, "m must be >= 1"
    assert 1 <= c <= p ** m, "column index out of range 1..p**m"
    val = c - 1
    x = []
    for _ in range(m):
        x.append(val % p)
        val //= p
    assert val == 0, "c-1 did not fit in m base-p digits"
    return tuple(x)


def point_to_column(x, p) -> int:
    """Inverse of column_to_point; 1-indexed column (NOTES Appendix C).

    c = 1 + sum_i x_i p^(i-1), with x_1 least significant; m = len(x).
    """
    assert p >= 2, "p must be >= 2"
    c_minus_1 = 0
    place = 1
    for xi in x:
        xi = int(xi)
        assert 0 <= xi < p, "coordinate out of range 0..p-1"
        c_minus_1 += xi * place
        place *= p
    return c_minus_1 + 1


def manhattan_weight(x) -> int:
    """Manhattan weight |x|_M = sum(int(x_i)) (NOTES sec 5, eq for W_M)."""
    return sum(int(xi) for xi in x)


def manhattan_puncture_columns(m, w, p) -> list:
    """Sorted 1-indexed columns with |x|_M <= w, i.e. the set S_p(w) (NOTES sec 5).

    Returns the column indices (1-based) of every point x in F_p^m whose
    Manhattan weight is at most w. The size equals [m,<=w]_p == pnomial_le(m,w,p).
    """
    assert p >= 2, "p must be >= 2"
    assert m >= 1, "m must be >= 1"
    cols = []
    for c in range(1, p ** m + 1):
        if manhattan_weight(column_to_point(c, m, p)) <= w:
            cols.append(c)
    return cols
