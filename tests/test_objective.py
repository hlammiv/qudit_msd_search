"""Protocol-cost search objective (leading-order surrogate) vs exact curves."""
import pytest

from qmsd.cyclic import construction_B
from qmsd.objective import (LeadingOrderCurve, marginal_saving, surrogate_from_matrices,
                            surrogate_from_structured_ad)
from qmsd.performance import leading_order, profile_exact, record_matrices
from qmsd.protocol import _rm_k1, curve
from qmsd.results import load_result_catalog


@pytest.fixture(autouse=True, scope="module")
def _tmp_hist_dir(tmp_path_factory):
    import qmsd.protocol as proto

    old = proto.HIST_DIR
    proto.HIST_DIR = tmp_path_factory.mktemp("hist")
    yield
    proto.HIST_DIR = old


@pytest.mark.parametrize("aid", ["paper-table3-p5-m2-n20-k5", "paper-table3-p3-m4-n72-k9"])
def test_structured_ad_gives_the_exact_per_output_count(aid):
    rec = {r.artifact_id: r for r in load_result_catalog()}[aid]
    sur = surrogate_from_structured_ad(rec.p, rec.m, rec.r, rec.puncture_columns)
    exact = leading_order(profile_exact(*record_matrices(rec), rec.p, kinds=("avg",)), "avg")
    assert (sur.d, sur.k) == (exact.d, rec.k)
    assert sur.L_d == pytest.approx(float(exact.L_d))


@pytest.mark.parametrize("aid", ["paper-table3-p5-m2-n20-k5", "imported-p7-n44-k5-1c67ab44c96c"])
def test_matrix_surrogate_matches_structured_ad(aid):
    """The generic (G0, H1) counter, needed for weighted supports, agrees with structured_ad."""
    rec = {r.artifact_id: r for r in load_result_catalog()}[aid]
    sur = surrogate_from_structured_ad(rec.p, rec.m, rec.r, rec.puncture_columns)
    for method in ("mitm", "balanced"):
        gen = surrogate_from_matrices(*record_matrices(rec), rec.p, method=method)
        assert (gen.n, gen.k, gen.d) == (sur.n, sur.k, sur.d)
        assert gen.L_d == pytest.approx(sur.L_d)


def test_surrogate_marginal_saving_tracks_exact():
    base = [curve("[[6,1,2]]_7 RS", *_rm_k1(7, 1), 7, key="obj-rs7")]
    exact = curve("[[7,1,3]]_7 cyc-B", *construction_B(7, 8, {1, 4, 7}), 7, key="obj-cyc7")
    sur = LeadingOrderCurve("surrogate", 7, exact.n, exact.k, exact.d, exact.c_eps * 6 ** exact.d)
    grid = ((1e-2, 1e-10), (1e-3, 1e-15))
    e, s = marginal_saving(exact, base, grid)["best"], marginal_saving(sur, base, grid)["best"]
    assert e > 1 and s == pytest.approx(e, rel=0.05)
