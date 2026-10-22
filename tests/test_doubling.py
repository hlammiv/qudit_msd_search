"""Qudit doubling and the Prakash-Saha d=2 codes, against arXiv:2610.29457 Tables 3 and 6."""
import numpy as np
import pytest

from qmsd.doubling import (PCODE_10_3, PCODE_7W_3, PCODE_7W_7, is_pcode, paper_codes, pcode_distance,
                           prakash_saha_d2)
from qmsd.performance import leading_order
from qmsd.protocol import best_protocol, curve
from qmsd.triorthogonal import is_triorthogonal_cubic


@pytest.fixture(autouse=True, scope="module")
def _tmp_hist_dir(tmp_path_factory):
    import qmsd.protocol as proto

    old = proto.HIST_DIR
    proto.HIST_DIR = tmp_path_factory.mktemp("hist")
    yield
    proto.HIST_DIR = old


def test_pcodes_of_table1():
    assert is_pcode(PCODE_7W_3[0], 3, PCODE_7W_3[1]) and pcode_distance(PCODE_7W_3[0], 3) == 3
    assert is_pcode(PCODE_10_3, 3) and pcode_distance(PCODE_10_3, 3) == 4
    assert is_pcode(PCODE_7W_7[0], 7, PCODE_7W_7[1]) and pcode_distance(PCODE_7W_7[0], 7) == 4
    assert not is_pcode(PCODE_7W_3[0], 3)          # the [[7,1,3]] code needs its charges


# label -> (d, A_d, c_eps, eps*, P_acc(0.01)) from Table 3
TABLE3 = {"[[15,1,3]]_3 dbl": (3, 50, 6.25, 0.230, 0.86),
          "[[28,1,4]]_3 dbl": (4, 230, 14.4, 0.212, 0.75),
          "[[37,1,5]]_3 dbl": (5, 2002, 62.6, 0.230, 0.69),
          "[[33,1,5]]_5 dbl": (5, 2552, 2.49, 0.334, 0.72),
          "[[20,1,4]]_7 dbl": (4, 1194, 0.921, 0.376, 0.82)}


@pytest.mark.parametrize("p", [3, 5, 7])
def test_doubled_codes_reproduce_table3(p):
    for label, G0, H1, coeffs in paper_codes(p):
        assert is_triorthogonal_cubic(np.vstack([H1, G0]), 1, p, coeffs)
        c = curve(label, G0, H1, p, key=f"test-{label}", coeffs=coeffs)
        d, A_d, c_eps, thr, pacc = TABLE3[label]
        assert c.d == d and leading_order(c.profile, "avg").L_d == A_d
        assert round(c.c_eps, 3 if c_eps < 1 else 2 if c_eps < 10 else 1) == c_eps
        assert round(float(c.threshold_eps), 3) == thr
        assert round(float(c.pacc(0.01)), 2) == pacc


def test_prakash_saha_matches_paper_example():
    G0, H1 = prakash_saha_d2(2)       # [[14,4,2]]_3, arXiv:2403.06228 eq. 11
    assert (H1 == np.array([[int(x) for x in r] for r in
                            ("11000000000222", "00110000000222", "00001100000222", "00000011000222")])).all()
    assert (G0 == np.array([[int(x) for x in r] for r in ("12121212012012", "00000000111222")])).all()
    for m in range(1, 13):
        G0, H1 = prakash_saha_d2(m)
        assert is_triorthogonal_cubic(np.vstack([H1, G0]), 3 * m - 2, 3)


def test_table6_costs_with_prakash_saha_codes_only():
    ps = [curve(f"PS m={m}", *prakash_saha_d2(m), 3, key=f"test-ps-{m}") for m in range(1, 13)]
    eps0 = 2 / 3 * 0.02                          # delta_in = 0.02, targets delta_out = 10^-x
    costs = [best_protocol(ps, eps0, 2 / 3 * 10.0 ** -x, max_depth=5).cost for x in (10, 15, 20, 25)]
    assert [round(c) for c in costs] == [49, 98, 110, 153]
