"""Tests for qmsd.asymptotics (NOTES sec 8, Appendix A; eqs 37/45/47)."""
import math

import pytest

from qmsd.asymptotics import H_p, gamma0, optimal_gamma, saddle_xi
from tests import ground_truth as gt


def _binary_entropy(theta):
    return -theta * math.log2(theta) - (1.0 - theta) * math.log2(1.0 - theta)


@pytest.mark.parametrize("theta", [0.05, 0.1, 0.25, 0.3, 0.4, 0.6, 0.75, 0.9])
def test_Hp_matches_binary_entropy(theta):
    # H_p is in log base p, so H_2 is the binary entropy in bits.
    assert H_p(theta, 2) == pytest.approx(_binary_entropy(theta), abs=1e-9)


def test_Hp_peak_is_one():
    # H_p peaks at theta = (p-1)/2 with value 1.
    for p in (2, 3, 5, 7, 11):
        assert H_p((p - 1) / 2.0, p) == pytest.approx(1.0, abs=1e-9)


def test_Hp_endpoints():
    for p in (2, 3, 5, 7):
        assert H_p(0.0, p) == pytest.approx(0.0, abs=1e-12)
        assert H_p(float(p - 1), p) == pytest.approx(0.0, abs=1e-12)


def test_saddle_p2_p3_closed_forms():
    # p=2 closed form: xi = theta/(1-theta); both satisfy theta(xi)=theta.
    for theta in (0.2, 0.45, 0.7):
        assert saddle_xi(theta, 2) == pytest.approx(theta / (1.0 - theta), rel=1e-9)
    # generic numeric saddle for p=5 must satisfy eq 45 residual ~ 0.
    from qmsd.asymptotics import _theta_of_xi

    for theta in (0.5, 1.5, 2.5, 3.5):
        xi = saddle_xi(theta, 5)
        assert _theta_of_xi(xi, 5) == pytest.approx(theta, abs=1e-8)


def test_gamma0_branch_continuity():
    # Both branches of eq 37 agree at the boundary theta = (p-1)/6.
    for p in (2, 3, 5, 7, 11):
        b = (p - 1) / 6.0
        left = gamma0(b - 1e-7, p)
        right = gamma0(b + 1e-7, p)
        assert left == pytest.approx(right, abs=1e-4)


def test_gamma0_requires_3theta_lt_pminus1():
    with pytest.raises(AssertionError):
        gamma0((2 - 1) / 3.0 + 0.01, 2)


@pytest.mark.parametrize("p", [2, 3, 5, 7, 11, 13, 17, 19, 23])
def test_optimal_gamma_reproduces_table1(p):
    g0, t0 = optimal_gamma(p)
    g_ref, t_ref = gt.TABLE1_ASYMPTOTIC[p]
    assert g0 == pytest.approx(g_ref, abs=2e-3)
    assert t0 == pytest.approx(t_ref, abs=2e-3)


def test_optimal_gamma_decreases_with_p():
    primes = sorted(gt.TABLE1_ASYMPTOTIC)
    vals = [optimal_gamma(p)[0] for p in primes]
    assert all(a > b for a, b in zip(vals, vals[1:]))
