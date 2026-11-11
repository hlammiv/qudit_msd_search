"""RigorousBoundCurve: f_lo <= surrogate <= f_hi, the bound dominates an exact curve where one
exists (RS code), the tail is the Singleton tail, and the optimizer's proved chain re-checks."""
import mpmath
import pytest

from qmsd.rigorous_curve import RigorousBoundCurve, from_census, rigorous_saving
from qmsd.objective import LeadingOrderCurve, standard_library
from qmsd.protocol import CodeCurve
from qmsd.rs_family import rs_profile, rs_parameters
from qmsd.performance import leading_order


def _rs_curve(p, k):
    prof = rs_profile(p, k); n, k, d = rs_parameters(p, k)
    return CodeCurve(f"[[{n},{k},{d}]]_{p} RS", p, n, k, prof)


def test_bounds_bracket_surrogate_and_exact_rs_curve():
    p, k = 13, 2
    ex = _rs_curve(p, k)
    L_d = ex.c_eps * (p - 1) ** ex.d
    # A_d of the RS magic code from its profile: use the leading order of the block (any A_d >= true
    # value keeps pacc_hi valid; take a generous one)
    rb = RigorousBoundCurve("rs", p, ex.n, ex.k, ex.d, L_d, A_d=10 ** 6)
    for eps in (1e-3, 1e-4, 1e-5, 1e-6):
        lo, hi, ex_f, sur = rb.f_lo(eps), rb.f_hi(eps), ex.f(eps), rb.f_surrogate(eps)
        assert lo <= sur * (1 + 1e-9) and sur <= hi
        assert ex_f <= hi, (eps, float(ex_f), float(hi))
        assert lo <= ex_f * (1 + 1e-6)


def test_tail_is_singleton_tail():
    rb = RigorousBoundCurve("x", 7, 50, 5, 3, L_d=1.0, A_d=6)
    eps = 1e-3
    e = mpmath.mpf(eps); eta = 7 * e / 6
    tail = sum(mpmath.binomial(50, w) * eta ** w * (1 - e) ** (50 - w) for w in range(4, 51)) * mpmath.mpf(7) ** (1 - 3)
    assert abs(rb._tail(eps) - tail) <= tail * 1e-12


def test_pacc_bounds_and_threshold():
    rb = RigorousBoundCurve("x", 7, 1756, 260, 7, L_d=1.29, A_d=24)
    for eps in (1e-4, 1e-3):
        assert rb.pacc_lo(eps) <= rb.pacc_hi(eps) <= 1.0
    thr = rb.threshold_eps
    assert thr is not None and 1e-3 < thr < 1e-2
    assert rb.f_hi(thr * 0.9) < thr * 0.9 and rb.f_hi(thr * 1.1) > thr * 1.1
    assert rb.exact is False


def test_from_census_and_proved_chain_rechecks():
    c = from_census("runs/weighted/p7_m4_sphere_1756_260_7_census.json")
    assert (c.n, c.k, c.d, c.A_d) == (1756, 260, 7, 24)
    res = rigorous_saving(c, standard_library(7), grid=((1e-2, 1e-30),))
    row = res["rows"][0]
    assert row["uses_proved"] and row["reaches_target_check"]
    assert row["saving_proved"] > 4 and row["saving_optimistic"] >= row["saving_proved"] - 1e-9


def test_correcting_curve_reduces_to_detection_and_trades_order_for_acceptance():
    from qmsd.rigorous_curve import CorrectingBoundCurve
    c0 = CorrectingBoundCurve("x", 11, 13332, 1309, 22, 0); c3 = CorrectingBoundCurve("x", 11, 13332, 1309, 22, 3)
    eps = 1e-4
    assert abs(c0.pacc_lo(eps) - (1 - eps) ** 13332) < 1e-9
    assert c3.pacc_lo(eps) > 0.9 > c0.pacc_lo(eps)
    assert c0.f_hi(eps) < c3.f_hi(eps) < 1e-12
    with pytest.raises(AssertionError):
        CorrectingBoundCurve("x", 7, 1745, 271, 7, 4)


def test_heavier_weights_reduce_to_singleton_and_tighten():
    """A heavier-weight bound equal to the Singleton count reproduces the plain curve; a smaller one tightens f_hi
    and pacc_hi (so f_lo can only rise)."""
    import mpmath
    from qmsd.rigorous_curve import RigorousBoundCurve, sphere_weight8_bound
    n, d, p = 1746, 7, 7
    plain = RigorousBoundCurve("c", p, n, 270, d, 193.2, 276)
    sing = int(mpmath.binomial(n, d + 1) * p ** 2)
    same = RigorousBoundCurve("c", p, n, 270, d, 193.2, 276, heavier={d + 1: (sing, sing)})
    for e in (1e-6, 1e-4, 3e-3):
        assert abs(float(same.f_hi(e)) / float(plain.f_hi(e)) - 1) < 1e-12
        assert abs(same.pacc_hi(e) - plain.pacc_hi(e)) < 1e-15
    N8 = sphere_weight8_bound(402)
    assert 9.4e16 < N8 < 9.5e16
    tight = RigorousBoundCurve("c", p, n, 270, d, 193.2, 276, heavier={d + 1: (N8, N8)})
    assert float(tight.f_hi(1e-6)) < 1e-3 * float(plain.f_hi(1e-6))
    assert float(tight.f_hi(1e-6)) > float(plain.f_lo(1e-6))
    for e in (1e-6, 3e-3):
        assert float(tight.f_lo(e)) >= float(plain.f_lo(e)) and tight.pacc_hi(e) <= plain.pacc_hi(e)
    assert tight.threshold_eps > plain.threshold_eps


def test_heavier_weights_validation():
    import pytest
    from qmsd.rigorous_curve import RigorousBoundCurve
    with pytest.raises(AssertionError):
        RigorousBoundCurve("c", 7, 100, 5, 7, 1.0, 6, heavier={9: (10, 10)})
    with pytest.raises(AssertionError):
        RigorousBoundCurve("c", 7, 100, 5, 7, 1.0, 6, heavier={8: (10, 11)})


def test_weight9_bound_tightens_further():
    from qmsd.rigorous_curve import RigorousBoundCurve, sphere_weight8_bound, sphere_weight9_bound, RM20_A9
    from math import comb
    N8, N9 = sphere_weight8_bound(402), sphere_weight9_bound(402)
    assert RM20_A9 == 9_202_538_394_000
    assert 1.0e20 < N9 < 1.2e20 and N9 < comb(1746, 9) * 7 ** 3 * 1e-5      # far below Singleton
    c8 = RigorousBoundCurve("c", 7, 1746, 270, 7, 193.2, 276, heavier={8: (N8, N8)})
    c9 = RigorousBoundCurve("c", 7, 1746, 270, 7, 193.2, 276, heavier={8: (N8, N8), 9: (N9, N9)})
    for e in (1e-6, 1e-5, 1e-4):
        assert float(c9.f_hi(e)) < float(c8.f_hi(e)) and float(c9.f_hi(e)) > float(c9.f_lo(e))
    assert c9.threshold_eps >= c8.threshold_eps
