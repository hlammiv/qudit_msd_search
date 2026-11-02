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
