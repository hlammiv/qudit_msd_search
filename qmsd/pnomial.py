"""p-nomial coefficients [m,s]_p and partial sums (NOTES sec 6).

The p-nomial coefficient is the Manhattan weight-distribution count

    (1 + x + x^2 + ... + x^{p-1})^m = sum_{s=0}^{m(p-1)} [m,s]_p x^s,
    [m,s]_p = #{ v in F_p^m : |v|_M = sum_i v_i = s }       (NOTES eqs 12/16).

Reduces to the ordinary binomial coefficient binom(m,s) when p=2.
Defined for integers m,s >= 0 and zero outside 0 <= s <= m(p-1).

Identities used here:
  - Lemma 3 (generalized Pascal, eq 19):  [m+1,s]_p = sum_{i=0}^{p-1} [m,s-i]_p.
  - Lemma 4 (multinomial, eq 20):         [m,s]_p = sum over compositions.
  - Lemma 5 (completeness, eq 22):        [m,>k]_p + [m,<=k]_p = p**m.
Cumulative counts (eqs 13/14/21):
  - [m,<=k]_p = number of punctures = code dimension k of weight-<=k puncturing,
  - [m,>k]_p  = surviving block length n.
"""
from __future__ import annotations

import math
from functools import lru_cache


@lru_cache(maxsize=None)
def _pnomial_row_cached(m, p) -> tuple:
    """Memoized core: coeffs of (1+x+...+x^{p-1})^m as a tuple of big ints."""
    assert p >= 2, "p must be >= 2"
    assert m >= 0, "m must be >= 0"
    row = (1,)                      # (...)^0 = 1
    base = (1,) * p
    for _ in range(m):
        row = _convolve(row, base)
    assert len(row) == m * (p - 1) + 1
    assert sum(row) == p ** m, "p-nomial row must sum to p**m"
    return row


def _convolve(a, b) -> tuple:
    """Discrete convolution of two integer sequences (polynomial multiply)."""
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return tuple(out)


def pnomial_row(m, p) -> list:
    """Coeffs of (1+x+...+x^{p-1})^m, length m(p-1)+1, by convolution; memoized.

    Repeatedly convolves the length-p all-ones vector m times (NOTES sec 6,
    eqs 12/16). Returns a fresh list of python ints (exact big-int arithmetic).
    """
    return list(_pnomial_row_cached(m, p))


def pnomial(m, s, p) -> int:
    """[m,s]_p; 0 if s<0 or s>m(p-1) (NOTES sec 6, eqs 12/16)."""
    if s < 0 or s > m * (p - 1):
        return 0
    return _pnomial_row_cached(m, p)[s]


def pnomial_le(m, s, p) -> int:
    """sum_{j<=s} [m,j]_p = [m,<=s]_p; 0 if s<0; p**m if s>=m(p-1).

    C_<=(m,s;W_M,p): number of punctures / code dimension of the weight-<=s
    Manhattan puncturing (NOTES eqs 13/21).
    """
    if s < 0:
        return 0
    top = m * (p - 1)
    if s >= top:
        return p ** m
    return sum(_pnomial_row_cached(m, p)[: s + 1])


def pnomial_gt(m, s, p) -> int:
    """p**m - pnomial_le(m,s,p) = [m,>s]_p.

    C_>(m,s;W_M,p): surviving block length n of the weight-<=s Manhattan
    puncturing (NOTES eqs 14/22, Lemma 5). For s<0 returns p**m.
    """
    return p ** m - pnomial_le(m, s, p)


def pnomial_multinomial(m, s, p) -> int:
    """[m,s]_p via Lemma 4 multinomial expansion (independent cross-check).

    Lemma 4 (NOTES eq 20):
        [m,s]_p = sum_{ k_1 + 2 k_2 + ... + (p-1) k_{p-1} = s,  k_i >= 0,
                        k_1 + ... + k_{p-1} <= m }
                    m! / ( k_1! ... k_{p-1}! (m - sum_i k_i)! ),
    where k_i counts coordinates equal to i (the rest equal 0). Computed by
    recursion over i = 1..p-1; used only to validate pnomial()/pnomial_row().
    """
    if s < 0 or s > m * (p - 1):
        return 0

    def recurse(i, remaining_sum, remaining_coords):
        if i == p:
            # Leftover coordinates are the zeros; weight must be exhausted.
            return 1 if remaining_sum == 0 else 0
        total = 0
        max_ki = min(remaining_coords, remaining_sum // i)
        for ki in range(max_ki + 1):
            sub = recurse(i + 1, remaining_sum - i * ki, remaining_coords - ki)
            if sub:
                total += math.comb(remaining_coords, ki) * sub
        return total

    return recurse(1, s, m)
