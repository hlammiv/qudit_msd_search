"""Decisive oracle test for the exact MacWilliams weight-distribution engine.

This pins ``qmsd.weightdist.exact_distance_and_Ad`` against the paper's published
ground truth (arXiv:2510.10852, Table 3) for the six *small-dual* oracle codes whose
shortened generator ``G0`` is low-dimensional enough to enumerate exactly:

    [[20,5,2]]_5   dim(G0)=1    A_d=760
    [[72,9,3]]_3   dim(G0)=6    A_d=648
    [[112,13,3]]_5 dim(G0)=7    A_d=512
    [[200,43,3]]_3 dim(G0)=8    A_d=1700
    [[206,37,4]]_3 dim(G0)=14   A_d=880     (3**14 = 4.78M words; near the budget)
    [[667,62,4]]_3 dim(G0)=16   A_d=3972    (3**16 = 43.0M words; SLOW)

For each we assert the engine reproduces (a) the certified minimum distance ``d``,
cross-checked against the independent meet-in-the-middle ``mindist.min_dependent_columns``;
and (b) the paper's ``A_d`` via ``B_d`` (the MacWilliams count of all weight-d dual
codewords).  Where it is cheap we also confirm ``B_d`` equals the slow-but-correct LOGICAL
reference ``qmsd.distance.A_d_logical_Z`` -- the empirical resolution of the critical
subtlety (no weight-d stabilizers, so all minimum-weight dual codewords are logical and
``A_d_logical == B_d``).

Everything is EXACT-INTEGER: every count compared here is a python int.
"""
from __future__ import annotations

import os
from math import comb

import pytest

from qmsd import weightdist as wd
from qmsd.distance import A_d_logical_Z
from qmsd.mindist import min_dependent_columns
from qmsd.oracle import load_oracle
from qmsd.triorthogonal import build_triorthogonal_code


# Published Table-3 minimum-weight logical counts for the six small-dual oracle codes.
_PAPER_AD = {
    "[[20,5,2]]_5": 760,
    "[[72,9,3]]_3": 648,
    "[[112,13,3]]_5": 512,
    "[[200,43,3]]_3": 1700,
    "[[206,37,4]]_3": 880,
    "[[667,62,4]]_3": 3972,
}

# (label, max_words, slow).  dim(G0) drives the message count q**dim:
#   dim<=14 enumerate in well under the default 5M budget and run fast;
#   [[667,62,4]] needs 3**16 = 43.0M words -> raise the budget and mark slow.
_FAST_CASES = [
    ("[[20,5,2]]_5", 5_000_000),
    ("[[72,9,3]]_3", 5_000_000),
    ("[[112,13,3]]_5", 5_000_000),
    ("[[200,43,3]]_3", 5_000_000),
    ("[[206,37,4]]_3", 5_000_000),
]
_SLOW_CASE = ("[[667,62,4]]_3", 50_000_000)

# A_d_logical_Z scans C(n,d) column subsets in pure python; only run it where that is cheap.
_SLOW_REF_SUBSET_CAP = 300_000


def _build_G0(label):
    oc = next(o for o in load_oracle() if o.label == label)
    built = build_triorthogonal_code(oc.p, oc.m, oc.r_max, oc.puncture_columns_1indexed)
    return oc, built, built["X_stab"]


def _check_oracle(label, max_words):
    oc, built, G0 = _build_G0(label)

    res = wd.exact_distance_and_Ad(G0, oc.p, max_words=max_words)

    # ---- certified minimum distance -------------------------------------------------
    assert res["feasible"] is True, label
    assert res["distance"] == oc.d, (label, res["distance"], oc.d)

    # Independent cross-check: distance of G0^perp = min # of F_p-dependent columns of G0.
    mdc = min_dependent_columns(G0, oc.p, d_max=oc.d)
    assert mdc == oc.d, (label, "mindist", mdc, oc.d)

    # ---- A_d (paper) reproduced exactly ---------------------------------------------
    paper_Ad = _PAPER_AD[label]
    assert res["B_d"] == paper_Ad, (label, "B_d", res["B_d"], paper_Ad)
    assert res["A_d_logical"] == res["B_d"], label
    assert res["A_d_equals_Bd"] is True, label

    # Sanity of the full dual weight distribution.
    B = res["weight_dist"]
    assert B[0] == 1 and B[oc.d] == paper_Ad
    assert all(B[w] == 0 for w in range(1, oc.d)), (label, "spurious sub-d weight")
    assert all(b >= 0 for b in B)

    # ---- empirical resolution of the logical-vs-all subtlety, where cheap -----------
    if comb(oc.n, oc.d) <= _SLOW_REF_SUBSET_CAP:
        ref = A_d_logical_Z(built, oc.p, oc.d)
        assert ref == res["B_d"] == paper_Ad, (label, ref, res["B_d"], paper_Ad)


@pytest.mark.parametrize("label,max_words", _FAST_CASES)
def test_oracle_distance_and_Ad_fast(label, max_words):
    _check_oracle(label, max_words)


@pytest.mark.slow
@pytest.mark.skipif(
    not os.environ.get("QMSD_RUN_SLOW"),
    reason="dim(G0)=16 -> 3**16 = 43.0M enumerated messages (minutes); set QMSD_RUN_SLOW=1",
)
def test_oracle_distance_and_Ad_667():
    label, max_words = _SLOW_CASE
    _check_oracle(label, max_words)


def test_infeasible_when_over_budget():
    # When q**dim(G0) exceeds max_words the engine must decline (no enumeration), not guess.
    _oc, _built, G0 = _build_G0("[[200,43,3]]_3")
    res = wd.exact_distance_and_Ad(G0, 3, max_words=10)
    assert res["feasible"] is False
    assert res["distance"] is None
    assert res["weight_dist"] is None
    assert res["B_d"] is None
    assert res["A_d_logical"] is None
    assert res["A_d_equals_Bd"] is None
