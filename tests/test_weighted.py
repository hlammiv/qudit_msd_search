"""Weighted Reed-Muller supports (charge distributions, arXiv:2610.29457 sec 4)."""
import numpy as np

from qmsd.triorthogonal import is_triorthogonal_cubic, logical_lambdas
from qmsd.weighted import (charge_degree_budget, cube_points, product_of_linear, weighted_code,
                           weighted_support)


def test_quadric_charges_give_a_shorter_triorthogonal_space():
    p, m = 7, 2                                   # r_max = 3, budget 12 - 1 - 9 = 2
    assert charge_degree_budget(p, m, 3) == 2
    g = product_of_linear(p, m, [((1, 0), 0), ((1, 0), 1)])     # x1 (x1 - 1)
    GX, X, mu = weighted_support(p, m, g)
    assert len(X) == 5 * 7
    K = GX.shape[0]
    assert is_triorthogonal_cubic(np.vstack([GX[:1], GX]), 0, p, mu)   # every triple vanishes
    # unit charges on the same support are NOT triorthogonal
    assert not is_triorthogonal_cubic(GX, 0, p)
    # puncture three cube-charge points -> a 3-output code with cube logical exponents
    pts = cube_points(mu, p)[:3]
    G0, H1, c = weighted_code(GX, mu, pts, p)
    assert is_triorthogonal_cubic(np.vstack([H1, G0]), 3, p, c)
    assert logical_lambdas(H1, 3, p, c) == [int(-mu[l] % p) for l in pts]
    assert H1.shape == (3, len(X) - 3) and G0.shape[0] == K - 3
