"""Tests for the structured A_d enumerator (qmsd.structured_ad).

Fast tests cover the F_p linear-algebra helpers, the flat enumeration count, and the
meet-in-the-middle codeword finder against brute force.  End-to-end tests cross-check the
structured A_d against the exact MacWilliams engine on the small-redundancy oracle codes,
and (guarded behind QMSD_SLOW=1, because they take minutes) reproduce the brute-force-
unreachable structured-only Table-3 values 572 / 1104 / 1128.
"""
import itertools
import os

import numpy as np
import pytest

from qmsd.oracle import load_oracle
from qmsd.reedmuller import r_max
from qmsd.structured_ad import (
    _flats,
    _mitm,
    _null_space,
    _rank,
    _solve,
    structured_ad,
)
from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.weightdist import exact_distance_and_Ad

SLOW = os.environ.get("QMSD_SLOW") == "1"
ORACLE = {c.label: c for c in load_oracle()}


def _gauss(n, k, q):
    num = den = 1
    for i in range(k):
        num *= q ** (n - i) - 1
        den *= q ** (k - i) - 1
    return num // den


# --------------------------------------------------------------------------- helpers
def test_null_space_and_solve_random():
    rng = np.random.default_rng(0)
    p = 3
    for _ in range(50):
        k, n = rng.integers(1, 6), rng.integers(2, 8)
        G = rng.integers(0, p, size=(k, n))
        N = _null_space(G, p)
        # every null vector v satisfies G v^T = 0
        assert N.shape[1] == n
        if N.shape[0]:
            assert np.all((G @ N.T) % p == 0)
        # rank-nullity
        assert _rank(G, p) + N.shape[0] == n


def test_solve_consistent():
    rng = np.random.default_rng(1)
    p = 3
    for _ in range(50):
        r, c = rng.integers(2, 8), rng.integers(1, 5)
        A = rng.integers(0, p, size=(r, c))
        u_true = rng.integers(0, p, size=c)
        b = (A @ u_true) % p
        u = _solve(A, b, p, c)
        assert np.all((A @ u) % p == b)


def test_flats_count():
    p = 3
    for m in (3, 4):
        for j in range(1, m + 1):
            flats = list(_flats(m, j, p))
            assert len(flats) == p ** (m - j) * _gauss(m, j, p)
            # each flat has p**j distinct points
            for cm in flats[:20]:
                assert len(set(int(x) for x in cm)) == p ** j


def _brute_weight_d_codewords(H, d, p):
    """All weight-d codewords of H^perp by brute column-subset scan (reference)."""
    rows, n = H.shape
    out = set()
    from qmsd.distance import _right_null_space, _full_support_kernel_vectors
    for T in itertools.combinations(range(n), d):
        sub = [H[:, j].tolist() for j in T]
        basis = _right_null_space(sub, p)
        for lam in _full_support_kernel_vectors(basis, p, d):
            out.add((T, tuple(lam)))
    return out


def test_mitm_matches_brute():
    rng = np.random.default_rng(2)
    p = 3
    for _ in range(30):
        rows, n = rng.integers(2, 5), rng.integers(5, 10)
        H = rng.integers(0, p, size=(rows, n))
        for d in (2, 3, 4):
            got = set(_mitm(H, d, p))
            exp = _brute_weight_d_codewords(H, d, p)
            assert got == exp, (rows, n, d)


# --------------------------------------------------------------------------- end-to-end
def test_structured_matches_macwilliams_72():
    """[[72,9,3]]_3 (m=4): full structured decomposition == exact MacWilliams A_d = 648.

    Fast (~6s) and exercises the whole machinery: j=2,3,4 flats, the b=1 (m=4) restricted
    code RM_3(1,2), the per-flat MITM, the lift and the minimal-span dedup.
    """
    c = ORACLE["[[72,9,3]]_3"]
    res = structured_ad(3, c.m, r_max(c.m, 3), c.puncture_columns_1indexed, jmax=c.m)
    assert res["distance"] == 3
    assert res["exact"] is True
    assert res["A_d"] == 648
    assert res["dim_histogram"] == {2: 504, 3: 144}     # no full-space (dim-4) tail
    # cross-check against the exact MacWilliams engine (dim G0 = 6, feasible)
    code = build_triorthogonal_code(3, c.m, r_max(c.m, 3), c.puncture_columns_1indexed)
    mac = exact_distance_and_Ad(code["G0"], 3, max_words=10 ** 8)
    assert mac["B_d"] == 648


def test_capped_jmax_is_certified_lower_bound():
    """A capped jmax never over-counts: structured A_d(jmax) <= true A_d, and is monotone."""
    c = ORACLE["[[72,9,3]]_3"]
    r = r_max(c.m, 3)
    a2 = structured_ad(3, c.m, r, c.puncture_columns_1indexed, jmax=2)
    a3 = structured_ad(3, c.m, r, c.puncture_columns_1indexed, jmax=3)
    assert a2["exact"] is False and a3["exact"] is False
    assert a2["A_d"] == 504                      # dim-2 only
    assert a2["A_d"] <= a3["A_d"] <= 648          # monotone, bounded by the true value


@pytest.mark.skipif(not SLOW, reason="set QMSD_SLOW=1 (minutes-long structured-only reproductions)")
def test_structured_only_230():
    """[[230,13,6]]_3: brute-force-unreachable (dim G0 = 38), structured A_d = 572."""
    c = ORACLE["[[230,13,6]]_3"]
    res = structured_ad(3, c.m, r_max(c.m, 3), c.puncture_columns_1indexed, jmax=3)
    assert res["distance"] == 6
    assert res["A_d"] == 572


@pytest.mark.skipif(not SLOW, reason="set QMSD_SLOW=1 (minutes-long structured-only reproductions)")
def test_structured_only_215():
    """[[215,28,5]]_3: brute-force-unreachable (dim G0 = 23), structured A_d = 1104."""
    c = ORACLE["[[215,28,5]]_3"]
    res = structured_ad(3, c.m, r_max(c.m, 3), c.puncture_columns_1indexed, jmax=4)
    assert res["distance"] == 5
    assert res["A_d"] == 1104


@pytest.mark.skipif(not SLOW, reason="set QMSD_SLOW=1 (minutes-long; needs the dim-5 tail)")
def test_brute_codes_exact_206():
    """[[206,37,4]]_3: full decomposition (jmax=m) reproduces the MacWilliams A_d = 880."""
    c = ORACLE["[[206,37,4]]_3"]
    res = structured_ad(3, c.m, r_max(c.m, 3), c.puncture_columns_1indexed, jmax=c.m)
    assert res["exact"] is True
    assert res["A_d"] == 880
