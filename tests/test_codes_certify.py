"""Integration tests for the exact MacWilliams path wired into qmsd.codes.

These pin the two capabilities the exact engine adds on top of the meet-in-the-middle
(MITM) minimum-distance routine, *through the codes.py API*:

  1. CAP LIFT (distance > 6).  ``mindist.min_dependent_columns`` certifies only d<=6.
     The exact engine certifies ANY distance from the dual weight distribution.  We use an
     MDS witness: a Reed-Solomon [n, kdim] generator over F_q has dual [n, n-kdim, kdim+1]
     (MDS), so its dual distance is kdim+1 -- choose kdim=6 -> distance 7, which the MITM
     provably cannot reach (it raises) but the engine certifies exactly.

  2. EXACT A_d for large codes.  ``A_d_logical_Z`` refuses when the C(n,d) subset scan
     exceeds its budget; the exact engine still returns A_d from the MacWilliams transform.

Everything is EXACT-INTEGER.
"""
from __future__ import annotations

import os

import numpy as np
import pytest

from qmsd import weightdist as wd
from qmsd.codes import code_certify, code_from_puncture
from qmsd.distance import A_d_logical_Z
from qmsd.mindist import min_dependent_columns
from qmsd.oracle import load_oracle
from qmsd.triorthogonal import build_triorthogonal_code


# ---------------------------------------------------------------------------
# 1. Cap lift: a code whose dual distance is 7 -- MITM cannot, the engine can.
# ---------------------------------------------------------------------------
def _reed_solomon_generator(q, n, kdim):
    """Vandermonde (Reed-Solomon) kdim x n generator over F_q (rows x^0..x^{kdim-1})."""
    pts = np.arange(n)
    return np.array(
        [[pow(int(x), i, q) for x in pts] for i in range(kdim)], dtype=np.int64
    ) % q


def test_macwilliams_certifies_distance_beyond_mitm_cap():
    # RS [10,6] over F_11 -> dual is MDS [10,4,7]: distance 7 > the MITM d<=6 cap.
    q, n, kdim = 11, 10, 6
    G0 = _reed_solomon_generator(q, n, kdim)

    res = wd.exact_distance_and_Ad(G0, q)
    assert res["feasible"] is True
    assert res["distance"] == 7, res["distance"]
    # MDS minimum-weight count A_{d} = C(n, d) (q-1):  C(10,7) * 10 = 120 * 10 = 1200.
    assert res["B_d"] == 1200
    assert res["weight_dist"][0] == 1
    assert all(res["weight_dist"][w] == 0 for w in range(1, 7))

    # The MITM cannot certify distance 7: it raises rather than mis-report.
    with pytest.raises(ValueError):
        min_dependent_columns(G0, q, d_max=8)


# ---------------------------------------------------------------------------
# 2. codes.py integration: code_certify reproduces paper d + A_d via the engine,
#    and the default code_from_puncture path is byte-for-byte unchanged.
# ---------------------------------------------------------------------------
_SMALL_DUAL = ["[[20,5,2]]_5", "[[72,9,3]]_3", "[[112,13,3]]_5", "[[200,43,3]]_3"]

_PAPER = {
    "[[20,5,2]]_5": (20, 5, 2, 760),
    "[[72,9,3]]_3": (72, 9, 3, 648),
    "[[112,13,3]]_5": (112, 13, 3, 512),
    "[[200,43,3]]_3": (200, 43, 3, 1700),
}


@pytest.mark.parametrize("label", _SMALL_DUAL)
def test_code_certify_matches_paper_and_default(label):
    oc = next(o for o in load_oracle() if o.label == label)
    n, k, d, A_d = _PAPER[label]

    cert = code_certify(oc.p, oc.m, oc.puncture_columns_1indexed)
    assert (cert.n, cert.k, cert.d, cert.A_d) == (n, k, d, A_d)
    assert cert.A_d_exact is True
    assert cert.d_certified is True

    # Default path (exact OFF) must agree -- no regression in the certified values.
    default = code_from_puncture(
        oc.p, oc.m, oc.puncture_columns_1indexed, max_distance=d + 1
    )
    assert (default.n, default.k, default.d, default.A_d) == (n, k, d, A_d)


def test_code_from_puncture_default_does_not_invoke_exact():
    # With exact_budget=0 (the default) the MacWilliams engine is never consulted; the
    # result is identical to the historical MITM/A_d_logical_Z path.  A_d_exact is still
    # set True here because A_d_logical_Z succeeded (small code) -- it marks an exact count,
    # not the engine specifically.
    oc = next(o for o in load_oracle() if o.label == "[[72,9,3]]_3")
    c = code_from_puncture(oc.p, oc.m, oc.puncture_columns_1indexed, max_distance=4)
    assert (c.n, c.k, c.d, c.A_d) == (72, 9, 3, 648)


# ---------------------------------------------------------------------------
# 3. Exact A_d where A_d_logical_Z refuses: [[667,62,4]] (C(667,4) ~ 8.2e9 > budget).
#    dim(G0)=16 -> 3**16 = 43.0M enumerated messages (minutes); gated behind QMSD_RUN_SLOW.
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not os.environ.get("QMSD_RUN_SLOW"),
    reason="3**16 = 43.0M enumerated messages (minutes); set QMSD_RUN_SLOW=1",
)
def test_code_certify_Ad_where_reference_refuses():
    oc = next(o for o in load_oracle() if o.label == "[[667,62,4]]_3")
    built = build_triorthogonal_code(oc.p, oc.m, oc.r_max, oc.puncture_columns_1indexed)

    # The slow logical reference must refuse this block (C(n,d) over budget).
    with pytest.raises(NotImplementedError):
        A_d_logical_Z(built, oc.p, oc.d)

    # The exact engine certifies the paper's A_d = 3972.
    cert = code_certify(
        oc.p, oc.m, oc.puncture_columns_1indexed, exact_budget=50_000_000
    )
    assert (cert.n, cert.k, cert.d, cert.A_d) == (667, 62, 4, 3972)
    assert cert.A_d_exact is True
