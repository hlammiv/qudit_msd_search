"""Tests for the memory-rebalanced weight-d counter (qmsd.weightcount_balanced).

The decisive correctness claim is the multiplicity normalization ``B_d = (p-1)*raw/C(d,d1)``
introduced by fixing the LEFT block's leading coefficient to 1.  These tests verify it by
EXACT agreement with TWO fully-independent references on small codes (several p, n, d), with
emphasis on the ``d = 6`` / ``d1 = 3`` path that the rebalance exists to make feasible:

  * ``count_weight_d``        -- the original full-coefficient MITM counter,
  * a brute column-subset count of all weight-d kernel words,
  * the exact ``MacWilliams`` dual weight distribution (computed from the code generator,
    independent of any MITM logic).

``count_weight_d_balanced(H, ...)`` counts weight-d codewords of ``ker(H)``; treating ``H``
as a generator, ``ker(H) = rowspace(H)^perp`` so ``macwilliams(weight_enumerator(H))[d]`` is
exactly the same quantity -- giving a third independent cross-check.
"""
import itertools
from math import comb

import numpy as np
import pytest

from qmsd.weightcount import count_weight_d
from qmsd.weightcount_balanced import count_weight_d_balanced, left_table_entries
from qmsd.weightdist import dual_weight_distribution


def _brute_count(H, d, p):
    """Number of weight-d codewords of ker(H) by brute column-subset scan (reference)."""
    from qmsd.distance import _right_null_space, _full_support_kernel_vectors

    rows, n = H.shape
    out = set()
    for T in itertools.combinations(range(n), d):
        sub = [H[:, j].tolist() for j in T]
        basis = _right_null_space(sub, p)
        for lam in _full_support_kernel_vectors(basis, p, d):
            out.add((T, tuple(lam)))
    return len(out)


@pytest.mark.parametrize("p", [3, 5, 7])
def test_balanced_matches_brute_and_count_weight_d(p):
    """balanced == count_weight_d == brute on random small parity checks (d up to 6)."""
    rng = np.random.default_rng(p)
    for _ in range(10):
        rows, n = int(rng.integers(2, 5)), int(rng.integers(7, 11))
        H = rng.integers(0, p, size=(rows, n))
        for d in (2, 3, 4, 5, 6):
            if d > n:
                continue
            bal = count_weight_d_balanced(H, d, p)
            assert bal == count_weight_d(H, d, p), (p, rows, n, d)
            assert bal == _brute_count(H, d, p), (p, rows, n, d)


# (p, n) with rows = n-3 -> ker dim ~3 (few codewords -> counters fast) AND p^rows small
# (MacWilliams feasible).  Covers the d=6 / d1=3 rebalanced path with B_6 > 0.
_MAC_CASES = [(3, 14), (3, 12), (5, 11), (5, 10), (7, 10), (7, 9), (11, 9), (11, 8), (13, 8)]


@pytest.mark.parametrize("p,n", _MAC_CASES)
def test_balanced_matches_macwilliams(p, n):
    """balanced == count_weight_d == MacWilliams B_d (independent of MITM), d up to 6."""
    rng = np.random.default_rng(1000 * p + n)
    rows = n - 3
    saw_d6_nonzero = False
    for _ in range(3):
        H = rng.integers(0, p, size=(rows, n))
        B = dual_weight_distribution(H, p, max_words=10 ** 8)  # B[d] = #weight-d of ker(H)
        for d in (2, 3, 4, 5, 6):
            bal = count_weight_d_balanced(H, d, p)
            assert bal == count_weight_d(H, d, p), (p, n, d, "vs count_weight_d")
            assert bal == B[d], (p, n, d, "vs MacWilliams", bal, B[d])
            if d == 6 and bal > 0:
                saw_d6_nonzero = True
    assert saw_d6_nonzero, f"p={p} n={n}: expected at least one B_6>0 case to exercise d1=3"


def test_p17_planted_weight6_nonzero():
    """At the pin's prime p=17, a PLANTED weight-6 dependency is counted exactly (B_6 > 0).

    Random columns are near-MDS (no weight-6 words), so we force cols {0..5} to sum to zero
    -> one weight-6 codeword and its 15 scalar multiples.  balanced == count_weight_d ==
    brute (brute is MITM-independent) at p=17, d=6 confirms the (p-1) multiplicity factor and
    the d1=3 fixed-leading-coefficient build at the actual pin prime.
    """
    p = 17
    saw_nonzero = False
    for seed in range(3):
        rng = np.random.default_rng(99000 + seed)
        H = rng.integers(1, p, size=(8, 12))
        H[:, 5] = (-(H[:, 0] + H[:, 1] + H[:, 2] + H[:, 3] + H[:, 4])) % p
        bal = count_weight_d_balanced(H, 6, p)
        assert bal == count_weight_d(H, 6, p), seed
        assert bal == _brute_count(H, 6, p), seed
        assert bal >= p - 1  # at least the planted orbit's p-1 scalar multiples
        saw_nonzero = True
    assert saw_nonzero


def test_left_table_is_p_minus_1_smaller():
    """The rebalanced d=6 LEFT table is exactly (p-1)x smaller than the full one."""
    p, n = 17, 237
    full = comb(n, 3) * (p - 1) ** 3          # count_weight_d's d1=3 table (OOM)
    bal = left_table_entries(n, 3, p)          # rebalanced (fits)
    assert bal == comb(n, 3) * (p - 1) ** 2
    assert full == bal * (p - 1)
    assert bal == 560_811_520


def test_rows_zero_full_space():
    """Empty parity check: ker = whole space, every weight-d pattern is a codeword."""
    H = np.zeros((0, 8), dtype=np.int64)
    for d in (2, 3, 6):
        assert count_weight_d_balanced(H, d, 5) == comb(8, d) * 4 ** d
        assert count_weight_d_balanced(H, d, 5) == count_weight_d(H, d, 5)


def test_no_codeword_returns_zero():
    """A full-rank H with no weight-d dependency returns 0 (e.g. identity-like, small d)."""
    H = np.eye(6, dtype=np.int64)  # columns independent: no nonzero kernel vector at all
    for d in (2, 3, 4):
        assert count_weight_d_balanced(H, d, 5) == 0
        assert count_weight_d_balanced(H, d, 5) == count_weight_d(H, d, 5)
