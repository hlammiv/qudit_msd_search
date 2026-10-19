"""Ring generalized triorthogonality: hosting, screens, induced logical — against the stored witnesses,
the Phase-1 p2_engine reference, dense leakage and the tri-block CCZ H3 screen."""

import random
import sys
from pathlib import Path

import numpy as np
import pytest

from primepower_msd import phase_tables as PT
from primepower_msd import rgt
from primepower_msd.ringlinalg import right_kernel_howell
from primepower_msd.symplectic import Z4Code, blocks, symp

ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / "primepower_msd"
CODE513 = [[1, 2, 1, 1, 2, 1, 0, 3, 3, 1], [1, 1, 1, 2, 1, 2, 1, 3, 0, 2],
           [3, 2, 0, 1, 2, 0, 1, 1, 1, 0], [2, 1, 2, 0, 2, 2, 1, 3, 1, 3]]


def _p2():
    sys.path.insert(0, str(ROOT / "scratchpad"))
    import p2_engine
    return p2_engine


def test_822_hosts_carrier_but_logical_does_not_escape():
    c = Z4Code.load(PKG / "witness_8_2_2_hosts_carrier.npy")
    L = rgt.TransversalLayer.uniform(8, np.array(PT.CARRIER_Q))
    assert rgt.hosts(c, L, "support").hosted and rgt.hosts(c, L, "cosets").hosted
    assert rgt.uniform_screens(c, PT.CARRIER_Q)["all"]
    P2 = _p2()
    rows = [list(map(int, r)) for r in c.rows]
    assert P2.hosts([r[:8] for r in rows if any(r[:8])], [r[8:] for r in rows if not any(r[:8])], 8)
    rep = rgt.induced_logical_report(c, L)
    assert rep["level"] == 3 and rep["escapes_A2"] is False and rep["odd_third_difference"] is False
    assert rep["table"] == [0, 4, 0, 4, 12, 4, 12, 4, 0, 12, 0, 12, 4, 4, 4, 4]
    par = rgt.axis_parity_report(c)
    assert par["odd_weights"] == [4, 4] and par["logical_escape_possible"] is False


def test_d3_codes_host_no_single_qudit_table():
    for code in (Z4Code.from_symplectic(CODE513, 5), Z4Code.load(PKG / "code_12_2_3_d3_k2.npy")):
        for q in PT.strict_level3_single(4):
            assert not rgt.uniform_screens(code, q)["S1_qweight"]
            assert not rgt.hosts(code, rgt.TransversalLayer.uniform(code.n, np.array(q)), "cosets").hosted


def test_phase5_witness_reproduced():
    sys.path.insert(0, str(ROOT / "scratchpad"))
    import cubic
    c = Z4Code.load(PKG / "code_10_4_3_gamma0834.npy")
    phi = np.load(ROOT / "scratchpad" / "witness_10_4_3_phi.npy")
    L = rgt.layer_from_poly_mod4(phi, cubic.monomials(10, 3), 10)
    assert rgt.hosts(c, L, "support").hosted
    v = rgt.hosts(c, L, "cosets")
    assert v.hosted and v.exact
    rep = rgt.induced_logical_report(c, L)
    assert rep["m"] == 4 and rep["level"] == 3
    assert rep["escapes_A2"] is False and rep["odd_third_difference"] is True
    assert L.degree_report()["protocol"] is False


def test_single_qudit_layer_matches_p2_engine_on_random_8_qudit_codes():
    """Exact Phase-1 M1 reference (brute support enumeration) vs rgt spanning-point hosting."""
    P2 = _p2()
    rng = random.Random(3)
    n = 8
    agree = 0
    for _ in range(60):
        # one X-generator with q-weight 0 (so S1 holds) plus random pure-Z rows orthogonal to it
        while True:
            a = [rng.randrange(4) for _ in range(n)]
            if any(a) and rgt.q_weight(PT.CARRIER_Q, a) == 0:
                break
        K = right_kernel_howell(np.array([a]), 4)
        zrows = [list(map(int, K[i])) for i in range(len(K)) if rng.random() < 0.5]
        code = Z4Code.from_css([a], zrows, n)
        L = rgt.TransversalLayer.uniform(n, np.array(PT.CARRIER_Q))
        ours = rgt.hosts(code, L, "support").hosted
        ref = P2.hosts([a], zrows, n)
        assert ours == ref
        agree += 1
    assert agree == 60


def _rand_iso(rng, n, r):
    rows, tries = [], 0
    while len(rows) < r and tries < 2000:
        tries += 1
        v = [rng.randrange(4) for _ in range(2 * n)]
        if any(v) and all(symp(v, g, n) == 0 for g in rows):
            rows.append(v)
    return rows


def test_cosets_mode_matches_dense_leakage_random():
    rng = random.Random(0)
    nrng = np.random.default_rng(0)
    strict = [np.array(q) for q in PT.strict_level3_single(4)]
    pairs = PT.sample_level3(2, 20, nrng)
    checked = 0
    for _ in range(60):
        n = rng.choice([3, 4])
        code = Z4Code.from_symplectic(_rand_iso(rng, n, rng.randrange(1, n + 1)), n)
        if not (code.free_quotient_ok and code.soft == 0):
            continue
        if rng.random() < 0.5:
            L = rgt.TransversalLayer.per_qudit([rng.choice(strict) if rng.random() < 0.7 else None for _ in range(n)])
        else:
            i, j = rng.sample(range(n), 2)
            L = rgt.TransversalLayer(n, 4, [None] * n, [(i, j, rng.choice(pairs))], [])
        assert rgt.hosts(code, L, "cosets").hosted == (rgt.diagonal_leakage(code, L) < 1e-8)
        checked += 1
    assert checked >= 20


def test_css_hosting_positive_cases_match_dense():
    rng = random.Random(5)
    strict = [np.array(q) for q in PT.strict_level3_single(4)]
    hosted = checked = 0
    for _ in range(80):
        n = rng.choice([3, 4])
        MX = [[rng.randrange(4) for _ in range(n)]]
        K = right_kernel_howell(np.array(MX), 4)
        MZ = [list(map(int, K[i])) for i in range(len(K)) if rng.random() < 0.6]
        try:
            code = Z4Code.from_css(MX, MZ, n)
        except AssertionError:
            continue
        if not (code.free_quotient_ok and code.soft == 0) or code.k_genuine == 0:
            continue
        L = rgt.TransversalLayer.per_qudit([rng.choice(strict) if rng.random() < 0.8 else None for _ in range(n)])
        h = rgt.hosts(code, L, "cosets").hosted
        assert h == rgt.hosts(code, L, "support").hosted == (rgt.diagonal_leakage(code, L) < 1e-8)
        hosted += h
        checked += 1
    assert checked >= 20


def test_triblock_ccz_hosting_equals_decisive_H1H3_and_dense():
    """Base code: X-stabilizers S, logical-X g, Z-stabilizers = (S + <g>)^perp (k0 = 1), three copies,
    across-block ring-CCZ.  Hosting must equal decisive.hosts_symbolic (H1-H3) and the dense leakage."""
    sys.path.insert(0, str(ROOT / "scratchpad"))
    import decisive
    rng = random.Random(1)
    done = dense = 0
    for _ in range(300):
        n0 = 2
        S = [[rng.randrange(4) for _ in range(n0)]]
        g = [rng.randrange(4) for _ in range(n0)]
        if not any(S[0]) or not any(g):
            continue
        MZ = right_kernel_howell(np.array([S[0], g]), 4)
        MZ = [list(map(int, z)) for z in MZ]
        try:
            base = Z4Code.from_css(S, MZ, n0)
        except AssertionError:
            continue
        tri = blocks(base, 3)
        L = rgt.TransversalLayer.across_block_ccz(n0)
        h = rgt.hosts(tri, L, "support").hosted
        assert h == decisive.hosts_symbolic(S, g, n0)
        if dense < 2:
            assert h == (rgt.diagonal_leakage(tri, L) < 1e-8)
            dense += 1
        done += 1
        if done >= 40:
            break
    assert done >= 20
    assert L.degree_report([0] * 2 + [1] * 2 + [2] * 2)["protocol"] is True
    assert L.degree_report()["protocol"] is False


def test_layer_structure_helpers():
    L = rgt.TransversalLayer.uniform(4, np.array(PT.CARRIER_Q))
    assert L.component_levels() == 3 and L.escapes_physically()
    L2 = rgt.TransversalLayer.cs_layer(4, [(0, 1, 1)])
    assert L2.component_levels() == 2 and not L2.escapes_physically()
    L3 = L.plus(L2)
    assert L3.phi([1, 1, 0, 0]) == (2 * PT.CARRIER_Q[1] + 4) % 16
