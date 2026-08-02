"""Phase-1 geometric (flat-occupancy) distance certifier: the j=2 min-weight upper bound is
EXACT for the paper's qutrit m=5 arc codes, and equals d_RM - max_2flat_occupancy(S) -- the
flat-occupancy law, now a certifier that needs no MITM weight tables (reaches d>=7)."""
import numpy as np
import pytest

from qmsd.oracle import load_oracle
from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.mindist import min_dependent_columns
from qmsd.structured_distance import (
    structured_distance,
    max_flat_occupancy,
    geometric_distance_upper,
    geometric_distance_dual,
    weight_hierarchy,
    flat_lower_bound,
)

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


def test_no_go_distance_at_most_6_for_all_k_ge_3():
    # RIGOROUS NO-GO: any 3 puncture points share a common 2-flat (they span <= 2 dims), whose
    # weight-9 indicator is a codeword, so |supp\S| <= 9 - 3 = 6. Hence d(S) <= 6 for EVERY k>=3
    # at qutrit m=5 -> [[230,13,6]] is distance-optimal and d>=7 is impossible. The two facts:
    #   (i) max_2flat_occupancy(S) >= 3 for every k>=3, and (ii) d_upper = d_RM - max_2flat <= 6.
    for S in [(1, 2, 3), (34, 61, 95, 140, 152), (5, 40, 88, 121, 199, 233)]:
        assert max_flat_occupancy(3, 5, S, 2) >= 3
        assert structured_distance(3, 5, 3, S)["d_upper"] <= 6


@pytest.mark.parametrize("oc", _QUTRIT_M5, ids=[oc.label for oc in _QUTRIT_M5])
def test_geometric_distance_dual_is_exact_on_qutrit_m5(oc):
    # The point-restricted dual-MITM geometric certifier (works where codeword enumeration OOMs,
    # e.g. p=7) reproduces the exact distance on the flat-binding qutrit m=5 oracle codes.
    assert geometric_distance_dual(oc.p, oc.m, oc.r_max, oc.puncture_columns_1indexed, j=2) == oc.d


@pytest.mark.parametrize("oc", _QUTRIT_M5, ids=[oc.label for oc in _QUTRIT_M5])
def test_geometric_distance_upper_is_tight_and_valid_on_qutrit_m5(oc):
    # qutrit m=5 is the tight regime (rtilde == (m-2)(p-1) == 6, so p^2 == d_RM == 9):
    # the screen value p^2 - max_2flat equals d_RM - max_2flat and is EXACT here.
    g = geometric_distance_upper(oc.p, oc.m, oc.r_max, oc.puncture_columns_1indexed)
    assert g == 9 - max_flat_occupancy(oc.p, oc.m, oc.puncture_columns_1indexed, 2)
    assert g >= oc.d                                   # a certified UPPER bound
    assert g == oc.d                                   # and tight in this regime


def test_geometric_distance_upper_returns_none_when_2flat_not_a_codeword():
    # When (m-2)(p-1) > rtilde the full-2-flat indicator is NOT a dual codeword, so there is no
    # valid 2-flat bound -> None (the search must then fall back to the MITM, never skip).
    # p=3, m=3, r=4: (m-2)(p-1)=2 > rtilde = m(p-1)-r-1 = 6-4-1 = 1.
    assert geometric_distance_upper(3, 3, 4, (1, 2, 5)) is None


def test_weight_hierarchy_second_weight_is_12():
    # Minimal-affine-span hierarchy of RM_3(6,5): span-2 (2-flat indicator) has weight 9 = d_RM;
    # the next span (span-3) has weight 12. The second weight w2=12 is the lower-bound input.
    assert weight_hierarchy(3, 5, 3) == {2: 9, 3: 12}


def test_flat_lower_bound_pins_distance_in_small_k_regime():
    # Certified lower bound min(d_RM - max_2flat, w2 - k) meets the Phase-1 upper bound in the
    # small-k regime -> EXACT distance with no MITM. Cross-checked here against the MITM.
    S = (34, 61, 95, 140, 152)  # k=5, max_2flat=3
    lo = flat_lower_bound(3, 5, 3, S, w2=12)      # min(9-3, 12-5) = min(6,7) = 6
    up = structured_distance(3, 5, 3, S)["d_upper"]
    G0 = np.asarray(build_triorthogonal_code(3, 5, 3, S)["X_stab"], dtype=int) % 3
    assert lo == up == min_dependent_columns(G0, 3, d_max=6) == 6
