"""Phase-1 geometric (flat-occupancy) distance certifier: the j=2 min-weight upper bound is
EXACT for the paper's qutrit m=5 arc codes, and equals d_RM - max_2flat_occupancy(S) -- the
flat-occupancy law, now a certifier that needs no MITM weight tables (reaches d>=7)."""
import numpy as np
import pytest

from qmsd.oracle import load_oracle
from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.mindist import min_dependent_columns
from qmsd.structured_distance import structured_distance, max_flat_occupancy

_QUTRIT_M5 = [oc for oc in load_oracle() if oc.p == 3 and oc.m == 5]


@pytest.mark.parametrize("oc", _QUTRIT_M5, ids=[oc.label for oc in _QUTRIT_M5])
def test_geometric_upper_bound_is_exact_and_matches_mitm(oc):
    res = structured_distance(oc.p, oc.m, oc.r_max, oc.puncture_columns_1indexed, jmax=2)
    # the j=2 term equals d_RM - max points on any 2-flat
    assert res["d_upper"] == res["d_RM"] - max_flat_occupancy(oc.p, oc.m, oc.puncture_columns_1indexed, 2)
    # and it is EXACT here: matches the paper distance AND the independent MITM
    G0 = np.asarray(build_triorthogonal_code(oc.p, oc.m, oc.r_max,
                                             oc.puncture_columns_1indexed)["X_stab"], dtype=int) % oc.p
    assert res["d_upper"] == oc.d == min_dependent_columns(G0, oc.p, d_max=6)


def test_geometric_is_an_upper_bound_never_under_reports():
    # d_upper exhibits an explicit flat codeword, so it is always >= the true distance
    # (a certified upper bound). Spot-check on the highest-distance oracle code.
    oc = next(o for o in _QUTRIT_M5 if o.label == "[[230,13,6]]_3")
    G0 = np.asarray(build_triorthogonal_code(oc.p, oc.m, oc.r_max,
                                             oc.puncture_columns_1indexed)["X_stab"], dtype=int) % oc.p
    d_true = min_dependent_columns(G0, oc.p, d_max=6)
    assert structured_distance(oc.p, oc.m, oc.r_max, oc.puncture_columns_1indexed)["d_upper"] >= d_true
