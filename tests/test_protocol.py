"""Multi-round protocol search over exact per-output curves."""
import mpmath
import pytest

from qmsd.protocol import CodeCurve, _rm_k1, best_protocol, curve
from qmsd.cyclic import construction_B


@pytest.fixture(autouse=True, scope="module")
def _tmp_hist_dir(tmp_path_factory):
    import qmsd.protocol as proto

    old = proto.HIST_DIR
    proto.HIST_DIR = tmp_path_factory.mktemp("hist")
    yield
    proto.HIST_DIR = old


@pytest.fixture(scope="module")
def p7():
    rs = curve("[[6,1,2]]_7 RS", *_rm_k1(7, 1), 7, key="test-rs7")
    cyc = curve("[[7,1,3]]_7 cyc-B", *construction_B(7, 8, {1, 4, 7}), 7, key="test-cyc7")
    return rs, cyc


def test_curve_matches_paper(p7):
    rs, cyc = p7
    assert (cyc.d, round(cyc.c_eps, 2), round(cyc.threshold_eps, 3)) == (3, 0.97, 0.446)
    assert (rs.d, round(rs.c_eps, 2), round(rs.threshold_eps, 3)) == (2, 2.5, 0.232)


def test_protocol_cost_is_product_of_rounds(p7):
    rs, cyc = p7
    best = best_protocol([rs], 1e-2, 1e-6, max_depth=4)
    assert all(c is rs for c in best.codes)
    cost, e = 1.0, mpmath.mpf(1e-2)
    for c in best.codes:
        cost *= c.n / (c.k * c.pacc(e))
        e = c.f(e)
    assert best.cost == pytest.approx(cost) and e <= 1e-6


def test_search_prefers_cheaper_code(p7):
    rs, cyc = p7
    both = best_protocol([rs, cyc], 1e-2, 1e-15, max_depth=6)
    only_rs = best_protocol([rs], 1e-2, 1e-15, max_depth=8)
    assert both.cost <= only_rs.cost


def test_non_cube_lambda_excluded():
    from qmsd.cyclic import construction_A

    assert curve("[[12,1,3]]_7", *construction_A(7, 12, {1, 2, 7}), 7, key="test-12-7") is None


def test_bounded_curve_brackets_exact():
    import numpy as np

    from qmsd.performance import record_matrices
    from qmsd.protocol import bounded_curve_from_record, protocol_bounds
    from qmsd.results import load_result_catalog

    rec = {r.artifact_id: r for r in load_result_catalog()}["paper-table3-p3-m4-n72-k9"]
    ex = curve(rec.label, *record_matrices(rec), 3, key="test-72-exact")
    bc = bounded_curve_from_record(rec, max_words=3 ** 3 * 9, pool=6)
    for e in np.logspace(-9, np.log10(0.4), 25):
        assert bc.f_lo(e) <= ex.f(e) * (1 + 1e-9) and ex.f(e) <= bc.f_hi(e) * (1 + 1e-9)
        assert bc.pacc_lo(e) <= ex.pacc(e) * (1 + 1e-9) <= bc.pacc_hi(e) * (1 + 1e-9) ** 2
    lo, hi = bc.threshold_bounds
    assert lo <= ex.threshold_eps and (hi is None or ex.threshold_eps <= hi)
    for eps0, tgt in [(1e-2, 1e-15), (1e-3, 1e-30)]:
        exact = best_protocol([ex], eps0, tgt).cost
        b = protocol_bounds([bc], eps0, tgt)
        assert b["lower"].cost <= exact * (1 + 1e-9) and exact <= b["upper"].cost * (1 + 1e-9)


def test_exact_curves_are_monotone(p7):
    import numpy as np

    for c in p7:
        es = np.logspace(-8, np.log10(1.5 * c.threshold_eps), 120)
        fs = [float(c.f(e)) for e in es]
        assert all(b >= a for a, b in zip(fs, fs[1:]))
