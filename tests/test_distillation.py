"""Tests for qmsd.distillation (NOTES eqs 38/39).

Anchored on the paper's worked single-round example, the ququint code
[[519,106,5]]_5 at delta_in = 1e-3 (NOTES sec 7):
    nbar_T = (0.9992)^519 * 106 ~ 70
    cost   = 519 / nbar_T       ~ 7.4
    delta_out = 2180/(70*4*625) * 1e-15 * (0.9992)^514 ~ 8e-18
"""
import math

import pytest

from qmsd import distillation as dist
from tests import ground_truth as gt

# [[519,106,5]]_5 with A_d = 2180  (TABLE3_CODES last entry).
DISTILLATION_EXAMPLE = next(
    c for c in gt.TABLE3_CODES
    if (c["p"], c["n"], c["k"], c["d"]) == (5, 519, 106, 5)
)
DELTA_IN = 1e-3


def test_example_present_in_ground_truth():
    c = DISTILLATION_EXAMPLE
    assert c["A_d"] == 2180


def test_nbar_T_exponent_is_n():
    """nbar_T raises survival to the power n (not n/k); NOTES sec 7 correction."""
    c = DISTILLATION_EXAMPLE
    n, k, p = c["n"], c["k"], 5
    survival = 1.0 - (p - 1) / p * DELTA_IN
    expected = survival ** n * k
    got = dist.nbar_T(n, k, p, DELTA_IN)
    assert got == pytest.approx(expected, rel=1e-12)
    # The n/k exponent would give a wildly different (much larger) value.
    wrong = survival ** (n / k) * k
    assert abs(got - wrong) > 1.0
    # Paper's worked figure ~70.
    assert got == pytest.approx(70.0, abs=2.0)


def test_cost_about_7_point_4():
    c = DISTILLATION_EXAMPLE
    n, k, p = c["n"], c["k"], 5
    nbar = dist.nbar_T(n, k, p, DELTA_IN)
    C = dist.cost(n, nbar)
    assert C == pytest.approx(7.4, abs=0.5)


def test_delta_out_avg_order_of_magnitude():
    c = DISTILLATION_EXAMPLE
    n, k, d, A_d, p = c["n"], c["k"], c["d"], c["A_d"], 5
    do = dist.delta_out_avg(n, k, d, A_d, p, DELTA_IN)
    # Within an order of magnitude of the paper's 8e-18.
    assert 8e-19 <= do <= 8e-17
    assert math.isclose(math.log10(do), math.log10(8e-18), abs_tol=1.0)


def test_delta_out_per_output_matches_eq38():
    p, d, A_d_i = 5, 5, 2180
    expected = A_d_i / ((p - 1) * p ** (d - 1)) * DELTA_IN ** d
    assert dist.delta_out_per_output(A_d_i, d, p, DELTA_IN) == pytest.approx(expected, rel=1e-12)


def test_cost_inverse_of_yield():
    assert dist.cost(519, 70.0) == pytest.approx(519 / 70.0, rel=1e-12)


def test_nbar_T_zero_noise_is_k():
    assert dist.nbar_T(519, 106, 5, 0.0) == 106


def test_delta_out_avg_consistent_with_per_output_form():
    """eq 39 = (A_d/nbar_T) * (eq-38 single-op value) * survival^(n-d)."""
    c = DISTILLATION_EXAMPLE
    n, k, d, A_d, p = c["n"], c["k"], c["d"], c["A_d"], 5
    nbar = dist.nbar_T(n, k, p, DELTA_IN)
    survival = 1.0 - (p - 1) / p * DELTA_IN
    base = dist.delta_out_per_output(A_d, d, p, DELTA_IN)  # A_d/((p-1)p^(d-1)) delta^d
    expected = base / nbar * survival ** (n - d)
    assert dist.delta_out_avg(n, k, d, A_d, p, DELTA_IN) == pytest.approx(expected, rel=1e-12)
