"""F_p arithmetic primitives (NOTES sec 2: field conventions; power-sum identity).

The single-variable power sum underlies Lemma 1 (RM dual) and Theorem 1
(triorthogonality): sum_{x in F_p} x^e vanishes for 0 <= e < p-1 and for any
e that is not a positive multiple of p-1 (NOTES eq surrounding line 117-119).
"""
from __future__ import annotations

from functools import lru_cache

import galois


def field_power_sum(p, a):
    """Return sum_{x in F_p} x^a  (mod p).

    By the power-sum identity (NOTES sec 2, around eq following Lemma 1):
    the sum equals -1 == p-1 (mod p) iff a > 0 and (p-1) | a, otherwise 0.

    Special case a == 0: x^0 == 1 for every x (using 0^0 == 1 by convention),
    so the sum is p ones == 0 (mod p).
    """
    assert isinstance(p, int) and p >= 2, "p must be an integer prime >= 2"
    assert isinstance(a, int) and a >= 0, "a must be a non-negative integer"
    if a > 0 and a % (p - 1) == 0:
        return p - 1
    return 0


@lru_cache(maxsize=None)
def GFp(p):
    """Return galois.GF(p), the F_p field class (cached thin wrapper) (NOTES sec 2).

    p must be prime for galois.GF to define a field; the assertion guards the
    global convention that all field arithmetic is mod a prime p.
    """
    assert isinstance(p, int) and p >= 2, "p must be an integer prime >= 2"
    return galois.GF(p)
