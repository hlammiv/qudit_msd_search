"""Unit tests for qmsd.puncture (NOTES Appendix C; sec 5)."""
import pytest

from qmsd.puncture import (
    column_to_point,
    point_to_column,
    manhattan_weight,
    manhattan_puncture_columns,
)
from qmsd.pnomial import pnomial_le
from qmsd.oracle import load_oracle


@pytest.mark.parametrize("p,m", [(2, 1), (2, 4), (3, 3), (5, 2), (3, 4)])
def test_column_point_round_trip_all_columns(p, m):
    seen = set()
    for c in range(1, p ** m + 1):
        x = column_to_point(c, m, p)
        assert len(x) == m
        assert all(0 <= xi < p for xi in x)
        assert point_to_column(x, p) == c
        seen.add(x)
    # Bijection: every point in F_p^m appears exactly once.
    assert len(seen) == p ** m


def test_hand_example_p5_m4():
    # c=2 -> c-1=1 = 1 + 0*p + 0*p^2 + 0*p^3 -> (1,0,0,0)
    assert column_to_point(2, 4, 5) == (1, 0, 0, 0)
    assert point_to_column((1, 0, 0, 0), 5) == 2
    # c=1 is always the origin.
    assert column_to_point(1, 4, 5) == (0, 0, 0, 0)
    # least-significant digit check: c=6 -> c-1=5 -> (0,1,0,0)
    assert column_to_point(6, 4, 5) == (0, 1, 0, 0)


def test_base_p_significance():
    # p=3, m=3: c-1 = 1*1 + 2*3 + 1*9 = 16 -> c=17 -> (1,2,1)
    assert column_to_point(17, 3, 3) == (1, 2, 1)
    assert point_to_column((1, 2, 1), 3) == 17


def test_manhattan_weight():
    assert manhattan_weight((0, 0, 0)) == 0
    assert manhattan_weight((1, 2, 1)) == 4
    assert manhattan_weight([4, 4]) == 8


@pytest.mark.parametrize("p,m", [(2, 4), (3, 3), (5, 2), (3, 4)])
def test_manhattan_set_sizes_match_pnomial(p, m):
    top = m * (p - 1)
    for w in range(-1, top + 2):
        cols = manhattan_puncture_columns(m, w, p)
        assert len(cols) == pnomial_le(m, w, p)
        # sorted, 1-indexed, within range
        assert cols == sorted(cols)
        assert all(1 <= c <= p ** m for c in cols)
        # membership is exactly the Manhattan ball
        for c in cols:
            assert manhattan_weight(column_to_point(c, m, p)) <= w


def test_manhattan_set_w0_is_origin():
    # only the origin point (column 1) has Manhattan weight 0
    assert manhattan_puncture_columns(4, 0, 3) == [1]


def test_oracle_columns_decode_in_range():
    # Every oracle puncture column must decode to a valid F_p^m point and
    # round-trip back, confirming the Appendix C convention is consistent.
    for oc in load_oracle():
        for c in oc.puncture_columns_1indexed:
            x = column_to_point(c, oc.m, oc.p)
            assert len(x) == oc.m
            assert point_to_column(x, oc.p) == c
