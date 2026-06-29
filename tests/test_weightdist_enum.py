"""Focused tests for qmsd.weightdist.weight_enumerator (chunked exact enumeration).

The weight enumerator is an EXACT-INTEGER object: every entry is a count, and a single
float would silently corrupt it.  These tests pin it against (a) an independent brute-force
enumeration, (b) the closed-form weight distribution of first-order Reed-Muller codes, and
(c) the documented over-budget guard.
"""
from __future__ import annotations

from itertools import product

import numpy as np
import pytest

from qmsd import weightdist as wd
from qmsd.field import GFp
from qmsd.reedmuller import rm_generator


# ---------------------------------------------------------------------------
# Independent brute-force reference (pure python ints)
# ---------------------------------------------------------------------------
def _brute_weight_enumerator(G, q):
    """Weight enumerator of the code spanned by G (r x n) by direct enumeration."""
    G = np.asarray(G, dtype=np.int64) % q
    r, n = G.shape
    hist = [0] * (n + 1)
    for msg in product(range(q), repeat=r):
        word = np.zeros(n, dtype=np.int64)
        for i, c in enumerate(msg):
            if c:
                word = (word + c * G[i]) % q
        hist[int(np.count_nonzero(word))] += 1
    return hist


def _random_full_rank_generator(rng, q, n, k):
    """Random full-rank k x n generator over F_q (k <= n)."""
    GF = GFp(q)
    while True:
        M = rng.integers(0, q, size=(k, n))
        if int(np.linalg.matrix_rank(GF(M % q))) == k:
            return M.astype(np.int64) % q


# ---------------------------------------------------------------------------
# weight_enumerator vs brute force on small codes
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "q,n,k",
    [(2, 7, 4), (2, 9, 5), (3, 6, 3), (3, 8, 4), (5, 5, 2), (5, 6, 3), (7, 4, 2)],
)
def test_weight_enumerator_matches_brute(q, n, k):
    rng = np.random.default_rng(2024 + q * 100 + n * 7 + k)
    for _ in range(5):
        G = _random_full_rank_generator(rng, q, n, k)
        got = wd.weight_enumerator(G, q)
        assert got == _brute_weight_enumerator(G, q)
        # Full-rank generator -> exactly q**k distinct codewords; the zero word has weight 0.
        assert sum(got) == q ** k
        assert got[0] == 1


def test_weight_enumerator_rank_deficient_uniform_multiple():
    # A rank-deficient generator counts each codeword with equal multiplicity q**(r-rank);
    # weight_enumerator enumerates messages, so it returns that uniform multiple (and brute
    # force, which also enumerates messages, agrees exactly).
    G = np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]], dtype=np.int64)  # rows sum to 0 mod 2
    got = wd.weight_enumerator(G, 2)
    assert got == _brute_weight_enumerator(G, 2)
    assert sum(got) == 2 ** 3  # all 8 messages enumerated, even though rank is 2


def test_weight_enumerator_returns_python_ints():
    rng = np.random.default_rng(11)
    G = _random_full_rank_generator(rng, 3, 7, 3)
    got = wd.weight_enumerator(G, 3)
    assert all(isinstance(v, int) and not isinstance(v, np.integer) for v in got)
    assert len(got) == G.shape[1] + 1


def test_weight_enumerator_empty_generator():
    # Zero rows -> the code is {0}: only the all-zero codeword of weight 0.
    G = np.zeros((0, 5), dtype=np.int64)
    got = wd.weight_enumerator(G, 3)
    assert got == [1, 0, 0, 0, 0, 0]


# ---------------------------------------------------------------------------
# Reproduces a known Reed-Muller weight distribution
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("m", [3, 4])
def test_weight_enumerator_first_order_rm(m):
    # First-order RM_2(1, m) is a [2^m, m+1, 2^(m-1)] code with closed-form weight enumerator
    # A_0 = 1, A_{2^m} = 1, A_{2^(m-1)} = 2^(m+1) - 2, all other A_w = 0.
    n = 2 ** m
    G = rm_generator(1, m, 2)
    assert np.asarray(G).shape == (m + 1, n)

    got = wd.weight_enumerator(G, 2)
    expected = [0] * (n + 1)
    expected[0] = 1
    expected[n] = 1
    expected[2 ** (m - 1)] = 2 ** (m + 1) - 2
    assert got == expected
    assert sum(got) == 2 ** (m + 1)


def test_weight_enumerator_rm_8_4_4_exact():
    # Explicit smallest case: RM_2(1,3) = [8,4,4], enumerator [1,0,0,0,14,0,0,0,1].
    G = rm_generator(1, 3, 2)
    assert wd.weight_enumerator(G, 2) == [1, 0, 0, 0, 14, 0, 0, 0, 1]


def test_weight_enumerator_ternary_rm():
    # Ternary RM_3(1, 2): a length-9 code over F_3; cross-check against brute force.
    G = rm_generator(1, 2, 3)
    n = np.asarray(G).shape[1]
    assert n == 9
    got = wd.weight_enumerator(G, 3)
    assert got == _brute_weight_enumerator(G, 3)


# ---------------------------------------------------------------------------
# over-budget guard
# ---------------------------------------------------------------------------
def test_weight_enumerator_raises_over_budget():
    rng = np.random.default_rng(99)
    G = _random_full_rank_generator(rng, 5, 6, 4)  # 5**4 = 625 messages
    with pytest.raises(ValueError):
        wd.weight_enumerator(G, 5, max_words=100)


def test_weight_enumerator_budget_boundary_ok():
    # q**r exactly at the budget must succeed (strictly-greater is the failure condition).
    rng = np.random.default_rng(100)
    G = _random_full_rank_generator(rng, 3, 6, 3)  # 3**3 = 27 messages
    got = wd.weight_enumerator(G, 3, max_words=27)
    assert sum(got) == 27
    with pytest.raises(ValueError):
        wd.weight_enumerator(G, 3, max_words=26)


def test_weight_enumerator_chunking_matches_single_pass():
    # Many small chunks vs. brute force -> identical exact histogram (no chunk-boundary loss).
    rng = np.random.default_rng(7)
    G = _random_full_rank_generator(rng, 3, 10, 6)  # 3**6 = 729 messages
    assert wd.weight_enumerator(G, 3) == _brute_weight_enumerator(G, 3)
