"""Scalar arithmetic over the cyclic ring Z_{2^k} (and general Z_n).

This is the ring layer the prime-only ``qmsd.field`` cannot serve: ``galois`` represents the
*field* GF(2^k), whereas a native d = 2^k qudit is the *ring* Z_{2^k} (a chain ring with zero
divisors and a cyclic additive group). The single load-bearing structural facts used downstream:

  * over Z_{2^k} (k >= 2) the squaring/cubing maps are NON-additive (unlike Frobenius x->x^2
    over GF(2^e)) — this is the algebraic lever that lets a single-qudit cubic/quadratic phase
    avoid the field collapse (see LITERATURE_MAP.md, GALOIS_QUDIT_INVESTIGATION.md);
  * the canonical lift Z_{2^k} -> {0,..,2^k-1} ⊂ Z carries a genuine integer value whose square
    lives mod 2^{k+2}, the precision a single-qudit level-3 phase needs.
"""

from __future__ import annotations


def val2(n: int) -> int:
    """2-adic valuation of n (number of factors of 2); val2(0) = +inf is returned as -1."""
    n = int(n)
    if n == 0:
        return -1
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def is_unit(a: int, n: int) -> bool:
    """True iff a is invertible in Z_n (gcd(a, n) == 1)."""
    from math import gcd

    return gcd(int(a) % int(n), int(n)) == 1


def units(n: int) -> list[int]:
    """The multiplicative group of units of Z_n."""
    return [a for a in range(n) if is_unit(a, n)]


def is_zero_divisor(a: int, n: int) -> bool:
    """True iff a is a nonzero zero-divisor in Z_n."""
    a = int(a) % int(n)
    return a != 0 and not is_unit(a, n)


def squaring_is_additive(d: int) -> bool:
    """Whether x -> x^2 is additive on Z_d, i.e. (a+b)^2 == a^2 + b^2 (mod d) for all a, b.

    True over GF(2^e) (Frobenius); FALSE over the ring Z_{2^k} for k >= 2 — the anti-collapse
    lever. Returns the boolean so the calibration step can assert non-additivity.
    """
    for a in range(d):
        for b in range(d):
            if ((a + b) ** 2) % d != (a * a + b * b) % d:
                return False
    return True


def cubing_is_additive(d: int) -> bool:
    """Whether x -> x^3 is additive on Z_d. FALSE over Z_{2^k} for k >= 2."""
    for a in range(d):
        for b in range(d):
            if ((a + b) ** 3) % d != (a ** 3 + b ** 3) % d:
                return False
    return True


def lift_square_mod(x: int, d: int, extra_bits: int = 2) -> int:
    """Integer square of the canonical lift of x in Z_d, reduced mod d * 2^extra_bits.

    For a single-qudit quadratic phase diag(zeta_{N}^{x^2}) the relevant precision is
    N = 2 * d * 2^{extra_bits-1}; this helper exposes the carry-carrying integer x^2 used to
    build that phase exactly.
    """
    xv = int(x) % int(d)
    return (xv * xv) % (int(d) * (2 ** int(extra_bits)))
