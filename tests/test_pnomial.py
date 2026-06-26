"""Unit tests for qmsd.pnomial (p-nomial coefficients, NOTES sec 6)."""
import math
import random

import pytest

from qmsd.pnomial import (
    pnomial,
    pnomial_row,
    pnomial_le,
    pnomial_gt,
    pnomial_multinomial,
)
from tests import ground_truth as gt

PRIMES = [2, 3, 5, 7, 11]


def test_reduces_to_binomial_for_p2():
    # [m,s]_2 == binom(m,s)  (NOTES eq 16: p-nomial reduces to binomial at p=2).
    for m in range(0, 12):
        for s in range(0, m + 2):
            assert pnomial(m, s, 2) == math.comb(m, s)


def test_row_sums_to_p_to_m():
    # sum_s [m,s]_p == p**m  (the generating function evaluated at x=1).
    for p in PRIMES:
        for m in range(0, 6):
            row = pnomial_row(m, p)
            assert len(row) == m * (p - 1) + 1
            assert sum(row) == p ** m
            assert sum(pnomial(m, s, p) for s in range(m * (p - 1) + 1)) == p ** m


def test_out_of_range_is_zero():
    for p in PRIMES:
        for m in range(0, 5):
            assert pnomial(m, -1, p) == 0
            assert pnomial(m, m * (p - 1) + 1, p) == 0


def test_multinomial_matches_convolution():
    # Lemma 4 (eq 20) computed independently must equal the convolution result.
    rng = random.Random(20510)
    for _ in range(300):
        p = rng.choice(PRIMES)
        m = rng.randint(0, 6)
        s = rng.randint(-1, m * (p - 1) + 1)
        assert pnomial_multinomial(m, s, p) == pnomial(m, s, p)


def test_completeness_lemma5():
    # Lemma 5 (eq 22): [m,>k]_p + [m,<=k]_p == p**m, for all k incl. edges.
    for p in PRIMES:
        for m in range(0, 5):
            for s in range(-2, m * (p - 1) + 3):
                assert pnomial_le(m, s, p) + pnomial_gt(m, s, p) == p ** m


def test_cumulative_edge_cases():
    for p in PRIMES:
        for m in range(0, 5):
            top = m * (p - 1)
            assert pnomial_le(m, -1, p) == 0
            assert pnomial_gt(m, -1, p) == p ** m
            assert pnomial_le(m, top, p) == p ** m
            assert pnomial_le(m, top + 5, p) == p ** m
            assert pnomial_gt(m, top, p) == 0
            # le is the running prefix sum of the row.
            assert pnomial_le(m, 0, p) == pnomial(m, 0, p) == 1


def test_generalized_pascal_lemma3():
    # Lemma 3 (eq 19): [m+1,s]_p == sum_{i=0}^{p-1} [m,s-i]_p.
    for p in PRIMES:
        for m in range(0, 6):
            for s in range(0, (m + 1) * (p - 1) + 1):
                lhs = pnomial(m + 1, s, p)
                rhs = sum(pnomial(m, s - i, p) for i in range(p))
                assert lhs == rhs


def test_symmetry():
    # The row is a palindrome: [m,s]_p == [m, m(p-1)-s]_p.
    for p in PRIMES:
        for m in range(0, 6):
            top = m * (p - 1)
            for s in range(0, top + 1):
                assert pnomial(m, s, p) == pnomial(m, top - s, p)


def test_cross_check_single_puncture_n():
    # Single-puncture codes have w=0: k=[m,<=0]_p=1 and n=[m,>0]_p=p**m-1,
    # matching gt.SINGLE_PUNCTURE_CODES (NOTES eqs 14/15).
    for c in gt.SINGLE_PUNCTURE_CODES:
        p, m, n, k = c["p"], c["m"], c["n"], c["k"]
        assert pnomial_le(m, 0, p) == k == 1
        assert pnomial_gt(m, 0, p) == n == p ** m - 1
