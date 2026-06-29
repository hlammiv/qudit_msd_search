"""Regression tests for M1.5 (Howell-form ring linear algebra) and the M2 machinery.

Run:  python -m pytest primepower_msd/tests/test_m1_5_m2.py -q
"""

import itertools

import numpy as np
import pytest

from primepower_msd.ringlinalg import howell_form, in_span, module_size, dual_module, right_kernel
from primepower_msd.ring_css import build_css, is_self_orthogonal, code_distance
from primepower_msd.ring_transversal import certify_distillation_code
from primepower_msd.hierarchy_search import diagonal_clifford, antidifference


def _brute_module(gens, n):
    G = [[int(x) % n for x in r] for r in gens] or [[0]]
    nc = len(G[0]); S = set()
    for co in itertools.product(range(n), repeat=len(G)):
        v = [0] * nc
        for c, row in zip(co, G):
            if c:
                for j in range(nc):
                    v[j] = (v[j] + c * row[j]) % n
        S.add(tuple(v))
    return S


# ----------------------------------------------------------------- M1.5: Howell form vs brute force
@pytest.mark.parametrize("n", [4, 8])
def test_howell_matches_brute_module(n):
    rng = np.random.default_rng(0)
    for _ in range(150):
        m, nc = rng.integers(1, 4), rng.integers(1, 5)
        M = rng.integers(0, n, size=(m, nc))
        brute = _brute_module(M, n)
        H = howell_form(M, n)
        hmod = _brute_module(H, n) if len(H) else {tuple([0] * nc)}
        assert brute == hmod
        assert module_size(M, n) == len(brute)
        for _ in range(4):
            v = tuple(rng.integers(0, n, size=nc))
            assert in_span(v, M, n) == (v in brute)


@pytest.mark.parametrize("n", [4, 8])
def test_dual_is_orthogonal_and_right_size(n):
    rng = np.random.default_rng(1)
    for _ in range(60):
        nc = int(rng.integers(2, 5))
        M = rng.integers(0, n, size=(int(rng.integers(1, 3)), nc))
        dual = dual_module(M, n)
        for y in dual:
            assert np.all((M @ np.array(y)) % n == 0)
        # |M| * |M^perp| == n^nc  (nondegenerate dot product over Z_n)
        assert module_size(M, n) * len(dual) == n ** nc


def test_known_howell_small():
    assert howell_form([[2, 4], [4, 0]], 8).tolist() == [[2, 4]]
    assert module_size([[2, 0], [0, 2]], 4) == 4


# ----------------------------------------------------------------- M2: the trivial transversal witness
def test_trivial_transversal_witness_d4():
    """The [[3,1]]_4 code with X-stabilizers (0,0,2),(0,2,0) admits a UNIFORM transversal level-3 gate
    (antidiff(S^2)) inducing a level-3 logical gate — but it is distance 1 (logical qudit unencoded)."""
    code = build_css([(0, 0, 2), (0, 2, 0)], 4)
    assert code is not None and code.k == 1 and code.cyclic
    assert code_distance(code) == 1                      # trivial: cannot distill
    U = antidifference(diagonal_clifford(4, 2, 0), 4)    # antidiff(S^2 Z^0), a level-3 gate
    cert = certify_distillation_code(code, U)
    assert cert["transversal"] and cert["logical_strict_level3"]


def test_self_orthogonality_check():
    assert is_self_orthogonal([(0, 0, 2), (0, 2, 0)], 4)
    assert not is_self_orthogonal([(1, 0, 0)], 4)        # 1*1 = 1 != 0 mod 4
