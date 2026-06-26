"""Unit tests for qmsd.reedmuller (NOTES sec 4; Thm 1, Thm 2; Lemma 2)."""
from __future__ import annotations

import numpy as np
import pytest

from qmsd import reedmuller as rm
from qmsd.field import GFp
from qmsd.oracle import load_oracle
from tests import ground_truth as gt


def test_r_max_matches_oracle():
    """r_max uses the m(p-1) form and agrees with OracleCode.r_max for every code."""
    for oc in load_oracle():
        assert rm.r_max(oc.m, oc.p) == oc.r_max == (oc.m * (oc.p - 1) - 1) // 3
        assert 3 * rm.r_max(oc.m, oc.p) < oc.m * (oc.p - 1)


@pytest.mark.parametrize("m,p", [(2, 3), (4, 3), (5, 3), (3, 5), (4, 5), (2, 11)])
def test_r_tilde_relation(m, p):
    r = rm.r_max(m, p)
    assert rm.r_tilde(r, m, p) == m * (p - 1) - r - 1


def test_d_rm_known_values():
    """Spot-check Theorem 2 integer formula on hand-computed cases."""
    assert rm.d_rm(0, 3, 3) == 27          # degree-0: whole space, distance p^m
    assert rm.d_rm(2, 2, 3) == 3
    assert rm.d_rm(5, 4, 3) == 6
    assert rm.d_rm(10, 4, 5) == 15
    assert rm.d_rm(2, 2, 2) == 1           # p=2 full degree -> distance 1


def test_d_rm_single_puncture_codes():
    """Lemma 2: single-puncture distance = d_RM(rtilde(r_max,m),m,p) - 1, end-to-end."""
    for c in gt.SINGLE_PUNCTURE_CODES:
        p, m, d = c["p"], c["m"], c["d"]
        r = rm.r_max(m, p)
        rt = rm.r_tilde(r, m, p)
        assert rm.d_rm(rt, m, p) - 1 == d, c


def test_points_base_p_digits():
    """points(m,p)[c] are the base-p digits of c with x_1 least significant."""
    p, m = 3, 3
    pts = rm.points(m, p)
    assert pts.shape == (p ** m, m)
    for c in range(p ** m):
        expected = [(c // (p ** j)) % p for j in range(m)]
        assert list(pts[c]) == expected


def test_monomial_exponents_count_and_bounds():
    """Count of exponent tuples equals number of generator rows; bounds respected."""
    p, m, r = 3, 4, 2
    exps = rm.monomial_exponents(m, r, p)
    # uniqueness + constraints
    assert len(set(exps)) == len(exps)
    for a in exps:
        assert len(a) == m
        assert all(0 <= ai <= p - 1 for ai in a)
        assert sum(a) <= r
    # deterministic / stable order
    assert exps == rm.monomial_exponents(m, r, p)


def test_rm_generator_shape_and_rank():
    """Generator has shape (#monomials, p**m) and is full row rank (basis)."""
    for (r, m, p) in [(2, 4, 3), (1, 2, 5), (5, 4, 3), (6, 2, 11)]:
        G = rm.rm_generator(r, m, p)
        exps = rm.monomial_exponents(m, r, p)
        assert G.shape == (len(exps), p ** m)
        # monomials are linearly independent => full row rank
        assert np.linalg.matrix_rank(G) == len(exps)


def test_rm_full_space_is_identity_dimension():
    """RM_p(m(p-1),m) spans the whole space F_p^{p^m} (dimension p^m)."""
    p, m = 3, 2
    G = rm.rm_generator(m * (p - 1), m, p)
    assert G.shape[0] == p ** m
    assert np.linalg.matrix_rank(G) == p ** m


def test_rm_generator_min_distance_matches_d_rm():
    """Exhaustive minimum nonzero weight of small RM codes equals d_rm (Thm 2)."""
    for (r, m, p) in [(1, 2, 3), (2, 2, 3), (1, 2, 5)]:
        G = rm.rm_generator(r, m, p)
        GF = GFp(p)
        k = G.shape[0]
        best = None
        # iterate all p^k codewords (small)
        for idx in range(1, p ** k):
            coeffs = []
            t = idx
            for _ in range(k):
                coeffs.append(t % p)
                t //= p
            cw = GF(np.array(coeffs, dtype=np.int64)) @ G
            w = int(np.count_nonzero(np.asarray(cw)))
            if w > 0 and (best is None or w < best):
                best = w
        assert best == rm.d_rm(r, m, p), (r, m, p)
