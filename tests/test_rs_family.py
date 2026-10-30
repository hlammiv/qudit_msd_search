"""Closed-form MDS performance of the Reed-Solomon family against the enumerating engine and the
published Table 2 values (arXiv:2610.08012)."""
import pytest

from qmsd.performance import canonical_logical_rows, leading_order, profile_exact, threshold
from qmsd.reedmuller import r_max, rm_generator
from qmsd.rs_family import mds_weight_enumerator, rs_profile, rs_parameters, r_max_rs
from qmsd.triorthogonal import build_triorthogonal_code


def _enumerated(p, k):
    r = r_max(1, p); G = rm_generator(r, 1, p)
    cols = list(range(1, k + 1))
    G0 = build_triorthogonal_code(p, 1, r, cols, G=G)["G0"]
    kinds = ("avg", "block") if p ** (r + 1) <= 5e7 else ("avg",)
    return profile_exact(G0, canonical_logical_rows(G, [c - 1 for c in cols], p), p, kinds=kinds)


@pytest.mark.parametrize("p,k", [(7, 1), (11, 1), (11, 2), (13, 3), (17, 2), (23, 5), (23, 4), (29, 7), (41, 10)])
def test_mds_profile_matches_enumeration(p, k):
    assert r_max_rs(p) == r_max(1, p)
    closed, enum = rs_profile(p, k), _enumerated(p, k)
    for eps in (0.3, 0.1, 0.02):
        assert float(closed.p_accept(eps)) == pytest.approx(float(enum.p_accept(eps)), rel=1e-12)
        for kind in enum.logical:
            assert float(closed.eps_out(eps, kind)) == pytest.approx(float(enum.eps_out(eps, kind)), rel=1e-9)
    lo_c, lo_e = leading_order(closed, "avg"), leading_order(enum, "avg")
    assert (lo_c.d, lo_c.L_d) == (lo_e.d, lo_e.L_d) and lo_c.d == rs_parameters(p, k)[2]


def test_table2_rs_codes():
    for p, c, thr, pacc in ((7, 2.5, 0.232, 0.94), (5, 1.5, 0.363, 0.96)):
        prof = rs_profile(p, 1)
        lo = leading_order(prof, "block")
        assert lo.c_eps == pytest.approx(c, abs=0.005)
        assert threshold(prof, "block").eps == pytest.approx(thr, abs=0.0015)
        assert float(prof.p_accept(0.01)) == pytest.approx(pacc, abs=0.005)


def test_mds_enumerator_sanity():
    A = mds_weight_enumerator(7, 2, 7)          # [7,2,6]_7
    assert A[6] == 7 * 6 and sum(A) == 49 and A[1:6] == [0] * 5
