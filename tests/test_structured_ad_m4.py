"""Tests for the m=4 / large-redundancy path of the structured A_d enumerator.

The structured enumerator's per-flat meet-in-the-middle (``_mitm``) used to encode each length-
``rows`` syndrome as one int64 p-adic, which raises ``OverflowError`` once ``p**rows > 2**63`` --
exactly the m=4, span-j>=3 regime (e.g. ``5**29``).  ``_mitm`` now hashes the syndrome (uint64,
overflow-proof) and verifies the true parity on each hash hit.  These tests pin:

  (1) the hashed ``_mitm`` matches brute force when ``p**rows`` far exceeds int64 (rows=30);
  (2) [SLOW] the full p=5 m=4 structured A_d on a high-puncture (dim G0 = 9) code equals the exact
      MacWilliams B_d -- a regime that is full-span (j=4) dominated and so genuinely exercises the
      large-redundancy span-4 flat that used to overflow.
"""
import itertools
import os

import numpy as np
import pytest

from qmsd.reedmuller import r_max, rm_generator
from qmsd.structured_ad import _mitm, _rank, structured_ad
from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.weightdist import exact_distance_and_Ad

SLOW = os.environ.get("QMSD_SLOW") == "1"


def _brute_weight_d_codewords(H, d, p):
    from qmsd.distance import _right_null_space, _full_support_kernel_vectors
    rows, n = H.shape
    out = set()
    for T in itertools.combinations(range(n), d):
        sub = [H[:, j].tolist() for j in T]
        basis = _right_null_space(sub, p)
        for lam in _full_support_kernel_vectors(basis, p, d):
            out.add((T, tuple(lam)))
    return out


def test_mitm_overflow_safe_large_rows():
    """_mitm is correct when p**rows >> int64 (rows=30, 5**30 ~ 9.3e20): hash encoder path."""
    p = 5
    rows, n = 30, 40
    assert p ** rows > np.iinfo(np.int64).max  # the old int64 p-adic would overflow here
    rng = np.random.default_rng(7)
    for _ in range(5):
        H = rng.integers(0, p, size=(rows, n))
        for d in (2, 3):
            got = set(_mitm(H, d, p))
            exp = _brute_weight_d_codewords(H, d, p)
            assert got == exp, (d,)


def _highk_fullrank_set(p, m, r, K, seed):
    """A size-K full-rank puncture set of RM_p(r,m) (so dim G0 = dim RM - K is small)."""
    Gi = np.asarray(rm_generator(r, m, p)).astype(np.int64) % p
    D = Gi.shape[0]
    rng = np.random.default_rng(seed)
    order = rng.permutation(p ** m)
    chosen, M = [], np.zeros((0, D), dtype=np.int64)
    for c in order:
        Mn = np.vstack([M, Gi[:, c]])
        if _rank(Mn, p) == Mn.shape[0]:
            M, _ = Mn, chosen.append(int(c))
        if len(chosen) == K:
            break
    return sorted(c + 1 for c in chosen)


@pytest.mark.skipif(not SLOW, reason="set QMSD_SLOW=1 (p=5 m=4 full structured run ~1 min)")
def test_structured_m4_p5_highk_matches_macwilliams():
    """p=5 m=4, K=113 punctures (dim G0 = 9): structured A_d (jmax=None) == MacWilliams B_d.

    These weight-3 dual codewords lift to FULL-SPAN (affine dim 4) RM_5(10,4) codewords, so the
    span-4 flat -- whose restricted-code redundancy is ``rows = 29`` (``5**29`` overflows int64) --
    is exactly what produces them; the count exercises the overflow-proof MITM end to end.
    """
    p, m = 5, 4
    r = r_max(m, p)
    cols = _highk_fullrank_set(p, m, r, 113, seed=113)
    built = build_triorthogonal_code(p, m, r, cols)
    G0 = np.asarray(built["G0"]) % p
    mac = exact_distance_and_Ad(G0, p, max_words=5_000_000)
    assert mac["feasible"] and mac["distance"] == 3
    res = structured_ad(p, m, r, cols, jmax=None, distance=3)
    assert res["A_d"] == mac["B_d"] == 680
    assert res["dim_histogram"] == {4: 680}  # all contributions are span-4 (full-space) codewords
