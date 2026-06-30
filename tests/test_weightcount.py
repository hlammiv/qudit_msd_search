"""Tests for the direct weight-d codeword counter (qmsd.weightcount.count_weight_d).

Cross-checks the streamed meet-in-the-middle count against the exact MacWilliams engine
(qmsd.weightdist.exact_distance_and_Ad) on the small-dual oracle codes, and against a brute
column-subset reference on random small parity checks.  The counter is dim(G0)-independent and
overflow-proof, so it reaches the regime where MacWilliams (p**dim G0) is infeasible.
"""
import itertools

import numpy as np
import pytest

from qmsd.oracle import load_oracle
from qmsd.reedmuller import r_max
from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.weightcount import count_weight_d
from qmsd.weightdist import exact_distance_and_Ad

ORACLE = {c.label: c for c in load_oracle()}


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


@pytest.mark.parametrize("p", [3, 5])
def test_count_matches_brute_random(p):
    rng = np.random.default_rng(p)
    for _ in range(15):
        rows, n = int(rng.integers(2, 5)), int(rng.integers(6, 11))
        H = rng.integers(0, p, size=(rows, n))
        for d in (2, 3, 4):
            if d > n:
                continue
            assert count_weight_d(H, d, p) == _brute_count(H, d, p), (p, rows, n, d)


@pytest.mark.parametrize("label", ["[[72,9,3]]_3", "[[206,37,4]]_3"])
def test_count_matches_macwilliams_oracle(label):
    """count_weight_d == MacWilliams B_d on oracle codes (the validated A_d)."""
    c = ORACLE[label]
    built = build_triorthogonal_code(3, c.m, r_max(c.m, 3), c.puncture_columns_1indexed)
    G0 = np.asarray(built["G0"]) % 3
    mac = exact_distance_and_Ad(G0, 3, max_words=10 ** 8)
    got = count_weight_d(G0, mac["distance"], 3)
    assert got == mac["B_d"]


def test_count_rows_zero_full_space():
    """Empty parity check: ker = whole space, every weight-d pattern is a codeword."""
    from math import comb
    H = np.zeros((0, 6), dtype=np.int64)
    assert count_weight_d(H, 2, 5) == comb(6, 2) * 4 ** 2
