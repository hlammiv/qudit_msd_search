"""Tests for the bounded m=3 structured-distance probe (qmsd/structured_m3.py).

Covers: AG(3,p) line/plane enumeration; validation 2a (the m=3 line logic restricted to a
plane reproduces structured_pe's m=2 lines); and validation 2b (d_struct == the certified
true distance) on fast p=5 / p=7 m=3 codes -- including a case where the plane (2-flat)
codeword caps the distance strictly below d_lines.  The slow p=11 m=3 checks (the actual
probe regime) are gated behind QMSD_SLOW=1.
"""
import os

import numpy as np
import pytest

from qmsd.reedmuller import r_max
from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.mindist import min_dependent_columns
from qmsd.weightdist import exact_distance_and_Ad
from qmsd.structured_pe import enumerate_lines as pe_lines
from qmsd.structured_m3 import (
    enumerate_lines_3d,
    enumerate_planes_3d,
    line_distance,
    plane_distance,
    flat_distance,
    flat_spread_puncture,
)
from qmsd.structured_ad import _rank

SLOW = os.environ.get("QMSD_SLOW") == "1"


# ---------------------------------------------------------------------------
# enumeration
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("p", [3, 5, 7, 11])
def test_line_count_and_affine_dim(p):
    lines = list(enumerate_lines_3d(p))
    assert len(lines) == (p * p + p + 1) * p * p
    # every line has p distinct points spanning affine dim 1; each point on p^2+p+1 lines
    from collections import Counter
    inc = Counter()
    for idx, _t in lines:
        assert len(set(idx)) == p
        P = np.array([[c % p, (c // p) % p, (c // p ** 2) % p] for c in idx])
        assert _rank((P[1:] - P[0]) % p, p) == 1
        for c in idx:
            inc[c] += 1
    assert set(inc.values()) == {p * p + p + 1}
    assert len(inc) == p ** 3


@pytest.mark.parametrize("p", [3, 5, 7, 11])
def test_plane_count_and_affine_dim(p):
    planes = list(enumerate_planes_3d(p))
    assert len(planes) == (p * p + p + 1) * p
    for idx, _cm in planes[: min(len(planes), 200)]:
        assert len(set(idx.tolist())) == p * p
        P = np.array([[c % p, (c // p) % p, (c // p ** 2) % p] for c in idx])
        assert _rank((P[1:] - P[0]) % p, p) == 2


# ---------------------------------------------------------------------------
# validation 2a: m=3 line logic restricted to a plane == structured_pe m=2 lines
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("p", [5, 7, 11])
def test_lines_in_plane_match_structured_pe(p):
    pe = set(frozenset(idx) for idx, _ in pe_lines(p))
    planes = list(enumerate_planes_3d(p))
    Fset = set(int(c) for c in planes[3][0].tolist())
    glob2int = {int(c): i for i, c in enumerate(planes[3][0].tolist())}
    m3 = set()
    for idx, _t in enumerate_lines_3d(p):
        if all(c in Fset for c in idx):
            m3.add(frozenset(glob2int[c] for c in idx))
    assert len(m3) == p * (p + 1)
    assert m3 == pe


# ---------------------------------------------------------------------------
# d_lines: closed form + explicit certificate is a genuine logical operator
# ---------------------------------------------------------------------------
def test_line_distance_certificate_in_kernel():
    p, m = 5, 3
    r = r_max(m, p)
    cols0, _ml = flat_spread_puncture(p, 12, seed=1, return_maxline=True)
    punc = [int(c) + 1 for c in cols0]
    code = build_triorthogonal_code(p, m, r, punc)
    X = np.asarray(code["X_stab"]).astype(np.int64) % p
    lin = line_distance(p, m, r, punc, X=X, n=code["n"])
    assert lin["beta_line"] == 0
    assert lin["d_lines"] == lin["d_RM"] - lin["max_line_punctures"]
    assert lin["certificate"]["in_kernel"] is True
    assert lin["certificate"]["weight"] == lin["d_lines"]


# ---------------------------------------------------------------------------
# validation 2b: d_struct == certified true distance (fast p=5/p=7 codes)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "p,k", [(5, 11), (5, 12), (5, 13), (7, 44), (7, 48)]
)
def test_validation_2b_fast(p, k):
    m = 3
    r = r_max(m, p)
    cols0 = flat_spread_puncture(p, k, seed=1)
    punc = [int(c) + 1 for c in cols0]
    code = build_triorthogonal_code(p, m, r, punc)
    assert code["full_rank"]
    G0 = np.asarray(code["G0"]).astype(np.int64) % p
    true_d = min_dependent_columns(G0, p, d_max=6)
    fd = flat_distance(p, m, r, punc, d_max=6)
    # d_struct (line + plane) must equal the certified true distance exactly
    assert fd["d_struct"] == true_d, (p, k, fd, true_d)
    # cross-check with MacWilliams when the dual is small enough to enumerate
    if p ** G0.shape[0] <= 5_000_000:
        mac = exact_distance_and_Ad(G0, p)["distance"]
        assert mac == true_d == fd["d_struct"]


def test_plane_can_cap_below_d_lines():
    """The decisive structure: a 2-flat codeword caps the distance strictly below d_lines."""
    p, m = 7, 3
    r = r_max(m, p)
    cols0 = flat_spread_puncture(p, 48, seed=1)
    punc = [int(c) + 1 for c in cols0]
    code = build_triorthogonal_code(p, m, r, punc)
    G0 = np.asarray(code["G0"]).astype(np.int64) % p
    true_d = min_dependent_columns(G0, p, d_max=6)
    lin = line_distance(p, m, r, punc)
    pl = plane_distance(p, m, r, punc, d_max=6, G0=G0,
                        S=set(int(c) for c in cols0))
    assert pl["d_planes"] == true_d
    assert pl["d_planes"] < lin["d_lines"]  # the plane (2D) codeword is the binding cap


# ---------------------------------------------------------------------------
# slow: p=11 m=3 (the actual probe regime) -- d_struct == certified true distance
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not SLOW, reason="set QMSD_SLOW=1 for the p=11 m=3 checks")
@pytest.mark.parametrize("k", [210, 208])
def test_validation_2b_p11(k):
    p, m = 11, 3
    r = r_max(m, p)
    cols0 = flat_spread_puncture(p, k, seed={210: 4, 208: 2}[k])
    punc = [int(c) + 1 for c in cols0]
    code = build_triorthogonal_code(p, m, r, punc)
    G0 = np.asarray(code["G0"]).astype(np.int64) % p
    true_d = min_dependent_columns(G0, p, d_max=6)
    # sound prune: a weight-<6 plane codeword is genuinely 2D (base weight >= 2(p-1)=20)
    pl = plane_distance(p, m, r, punc, d_max=true_d, G0=G0,
                        S=set(int(c) for c in cols0), only_heavy=2 * (p - 1) - true_d)
    lin = line_distance(p, m, r, punc)
    d_struct = min(lin["d_lines"], pl["d_planes"] if pl["d_planes"] else 99)
    assert d_struct == true_d
