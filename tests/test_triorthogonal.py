"""Unit tests for qmsd.triorthogonal (NOTES sec 3, Def 1, Thm 1)."""
from __future__ import annotations

import numpy as np
import pytest

from qmsd.field import GFp
from qmsd.reedmuller import r_max, rm_generator
from qmsd.triorthogonal import (
    build_triorthogonal_code,
    dual_matrix,
    is_triorthogonal,
    puncture_matrix,
    shorten_matrix,
)
from qmsd.oracle import load_oracle


# Small (p, m) pairs whose RM_p(r_max, m) generators are cheap to verify.
SMALL_PM = [(3, 2), (3, 3), (3, 4), (5, 2), (7, 2)]


@pytest.mark.parametrize("p, m", SMALL_PM)
def test_rm_at_rmax_is_triorthogonal(p, m):
    """RM_p(r_max, m) is triorthogonal since 3*r_max < m(p-1) (Thm 1)."""
    r = r_max(m, p)
    G = rm_generator(r, m, p)
    assert is_triorthogonal(G, p) is True


@pytest.mark.parametrize("p, m", SMALL_PM)
def test_over_degree_rm_not_triorthogonal(p, m):
    """RM_p(r_max+1, m) violates 3r < m(p-1), so it must fail triorthogonality."""
    r = r_max(m, p) + 1
    assert 3 * r >= m * (p - 1)          # genuinely over the threshold
    assert r <= m * (p - 1)              # still a valid RM degree
    G = rm_generator(r, m, p)
    assert is_triorthogonal(G, p) is False


def test_puncture_and_shorten_shapes_and_nesting():
    """G0 (shortened) is contained in Gp (punctured), and dual(Gp) is its X-dual."""
    p, m = 3, 3
    r = r_max(m, p)
    G = rm_generator(r, m, p)
    cols0 = [0, 5, 7]
    Gp = puncture_matrix(G, cols0, p)
    G0 = shorten_matrix(G, cols0, p)
    n = p ** m - len(cols0)
    assert Gp.shape[1] == n
    assert G0.shape[1] == n
    # G0 subset Gp: every shortened row lies in the rowspace of Gp.
    GF = GFp(p)
    stacked = GF(np.vstack([np.asarray(Gp), np.asarray(G0)]))
    assert int(np.linalg.matrix_rank(stacked)) == int(np.linalg.matrix_rank(Gp))
    # Dual is orthogonal to Gp.
    Z = dual_matrix(Gp, p)
    if Z.shape[0]:
        assert np.all(np.asarray(Gp @ Z.T) % p == 0)


def test_build_no_puncture_full_rank():
    """Empty puncture set: full rank, n = p**m, k = 0."""
    p, m = 3, 2
    code = build_triorthogonal_code(p, m, r_max(m, p), [])
    assert code["full_rank"] is True
    assert code["k"] == 0
    assert code["n"] == p ** m


@pytest.mark.parametrize("oc", load_oracle(), ids=lambda oc: oc.label)
def test_oracle_codes_build_full_rank(oc):
    """Every oracle code rebuilds to a full-rank [[n, k, d]] with n=p**m-k, k=#punctures."""
    code = build_triorthogonal_code(
        oc.p, oc.m, oc.r_max, oc.puncture_columns_1indexed
    )
    assert code["k"] == oc.num_punctures == len(oc.puncture_columns_1indexed)
    assert code["n"] == oc.p ** oc.m - code["k"]
    assert code["n"] == oc.n
    assert code["k"] == oc.k
    assert code["full_rank"] is True
