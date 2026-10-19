"""Lifts and the Galois-ring engine; anchors the first universality-enabling transversal ring codes."""

import numpy as np

from primepower_msd import lifts, rgt
from primepower_msd import phase_tables as PT
from primepower_msd.galoisring import gr, random_C0, random_symmetric_matrix
from primepower_msd.symplectic import Z4Code
from primepower_msd.z4_distance import code_distance


def test_gr_engine_reproduces_gr42_and_scales():
    R2 = gr(2)
    assert R2.h == [1, 1, 1] and R2.gram.tolist() == [[2, 3], [3, 3]] and R2.gram_det_unit
    for m in (2, 3, 4):
        R = gr(m)
        assert R.gram_det_unit
        assert len(R.teichmueller()) == 2 ** m
        assert R.pow(R.xi, 2 ** m) == R.xi                  # Frobenius^m = id on the Teichmueller root
        import random
        rng = random.Random(1)
        C0 = random_C0(R, 3, 2, rng); M = random_symmetric_matrix(R, 3, rng)
        rows, N = R.graph_code(C0, M)
        assert Z4Code.from_symplectic(rows, N).isotropic


def test_lifted_punctured_rm_codes():
    """[[15,1,3]]_4 hosts only non-escaping tables (X-stabilizer weight 8 != 0 mod 16);
    [[31,1,3]]_4 hosts all 48 and 32 induce an ESCAPING logical (weights 16): the first
    transversal universality-enabling ring-MSD code."""
    strict = [np.array(q) for q in PT.strict_level3_single(4)]
    expect = {4: (15, 16, 0), 5: (31, 48, 32)}
    for m, (n, n_hosted, n_escaping) in expect.items():
        b = lifts.punctured_rm_triorthogonal(m)
        c = lifts.lift_css_from_x(b["MX"], np.ones((1, n), dtype=np.int64), n)
        assert (c.n, c.k, c.k_genuine, c.soft, c.free_quotient_ok) == (n, 1.0, 1, 0, True)
        dmin, cert, det = code_distance(c, d_max=4)
        assert (dmin, cert, det["d_Z"]) == (3, True, 3)
        hosted = [q for q in strict if rgt.hosts(c, rgt.TransversalLayer.uniform(n, q), "support").hosted]
        assert len(hosted) == n_hosted
        esc = [q for q in hosted if PT.classify(rgt.induced_logical(c, rgt.TransversalLayer.uniform(n, q)))["escapes_A2"]]
        assert len(esc) == n_escaping
    # the carrier on [[31,1,3]] induces the strict level-3 escaping single-ququart table [0,15,12,3]
    b = lifts.punctured_rm_triorthogonal(5)
    c = lifts.lift_css_from_x(b["MX"], np.ones((1, 31), dtype=np.int64), 31)
    F = rgt.induced_logical(c, rgt.TransversalLayer.uniform(31, np.array(PT.CARRIER_Q)))
    assert F.tolist() == [0, 15, 12, 3] and PT.level(F) == 3


def test_cube_code_lift():
    cc = lifts.cube_8_3_2()
    c = lifts.lift_css(cc["MX"], cc["MZ"], 8)
    assert (c.k, c.k_genuine, c.soft) == (3.0, 3, 0)
    dmin, cert, det = code_distance(c, d_max=4)
    assert dmin == 2 and det == {"d_X": 4, "d_Z": 2}
    assert rgt.axis_parity_report(c)["logical_escape_possible"] is False


def test_random_generalized_triorthogonal_and_refinements():
    rng = np.random.default_rng(0)
    tri = lifts.random_binary_generalized_triorthogonal(16, 2, rng, n_stab=2, tries=60000)
    if tri is None:
        import pytest
        pytest.skip("no random generalized-triorthogonal matrix found in budget")
    G = np.vstack([tri["G1"], tri["G0"]])
    for a in range(4):
        for b in range(a, 4):
            if a == b and a < 2:
                assert G[a].sum() % 2 == 1
            else:
                assert np.sum(G[a] & G[b]) % 2 == 0
    c = lifts.lift_css_from_x(tri["G0"], tri["G1"], 16)
    assert c.isotropic and c.k_genuine <= 2
    import random
    refs = lifts.z_refinements(c, random.Random(0), n_samples=5)
    assert all(r.isotropic for r in refs)


def test_zrm_sizes_and_residue_bound():
    from primepower_msd.ringlinalg import module_size
    X = lifts.zrm(1, 3)
    # |ZRM(1,3)| = 4^{dim RM(0,3)} * 2^{dim RM(1,3) - dim RM(0,3)} = 4 * 2^3
    assert module_size(X, 4) == 4 * 8
    b = lifts.punctured_rm_triorthogonal(4)
    c = lifts.lift_css_from_x(b["MX"], np.ones((1, 15), dtype=np.int64), 15)
    rb = lifts.residue_distance_bound(c, d_max=6)
    assert rb["d_res_ub"] is not None and rb["d_res_ub"] >= 3
