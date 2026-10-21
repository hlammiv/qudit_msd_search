import numpy as np
import pytest

from primepower_msd.meas_distill import signed_rm
from primepower_msd.meas_sim import (Protocol, monte_carlo, leading_order, _closed_sets, closed_sets_fast,
                                     count_input_failures, count_input_failures_fast)

TOY_GRAPH = np.array([(a, b) for a in range(10) for b in range(a + 1, 10) if (b - a) % 10 not in (1, 9)])


def toy1(d=4):
    S, tau = signed_rm(2, 6)
    return Protocol(S, tau, d, np.repeat(np.arange(3), 21).reshape(63, 1))


def toy2(d=4):
    S, tau = signed_rm(1, 4)
    return Protocol(S, tau, d, TOY_GRAPH)


def test_toy1_pairs_count():
    # two cancelling inputs in one check: 3 checks * C(21,2) pairs * (d-1) value pairs
    L = leading_order(toy1(), 2)
    assert L["C"] == 3 * 210 * 3 and L["inner"] == 0


def test_toy2_weight3_matches_bruteforce():
    # 952 = numba brute force over all <= 3-fault patterns (scratch brute3.py; 0 below weight 3)
    P = toy2()
    assert leading_order(P, 3)["C"] == 952
    assert leading_order(P, 2)["C"] == 0


@pytest.mark.parametrize("W", [3, 4, 5])
def test_fast_enumeration_agrees(W):
    P = toy2()
    assert {x for x in _closed_sets(P, W) if len(x) >= 2} == closed_sets_fast(P, W)
    assert count_input_failures(P, W) == count_input_failures_fast(P, W)


def test_mc_acceptance_and_leading_order():
    P = toy1()
    p = 1e-3
    r = monte_carlo(P, p, shots=400_000, seed=11)
    acc = r["accepted"] / r["shots"]
    assert abs(acc - np.exp(-P.locations() * p)) < 0.01
    pred = 1890 * (p / 3) ** 2 * (1 - p) ** (P.locations() - 2) * r["shots"]
    assert 0.6 * pred < r["fail_shots"] < 1.5 * pred
