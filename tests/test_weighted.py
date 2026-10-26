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


def test_quadric_charges_spread_the_output_points():
    """p=7, m=3: the parabolic quadric keeps RM(5,3) triorthogonal on its 294-point support and
    its 98 cube-charge points meet every plane in at most 16 points (x1(x1-1): 49)."""
    import itertools
    from qmsd.reedmuller import points
    from qmsd.weighted import class_points, cube_points, quadric, weighted_support
    from qmsd.triorthogonal import is_triorthogonal_cubic
    p, m = 7, 3
    P = np.asarray(points(m, p), dtype=np.int64)
    GX, X, mu = weighted_support(p, m, quadric(p, m, (1, 1, 0), (0, 0, 1)))
    assert len(X) == 294 and is_triorthogonal_cubic(GX, 0, p, coeffs=mu)
    assert class_points(mu, p, 1) == cube_points(mu, p)
    def max_on_plane(X, idx):
        pts = P[X[idx]]
        return max(int(np.sum((pts @ np.asarray(a)) % p == b))
                   for a in itertools.product(range(p), repeat=3) if any(a) for b in range(p))
    assert len(cube_points(mu, p)) == 98 and max_on_plane(X, cube_points(mu, p)) == 16
    _, X2, mu2 = weighted_support(p, m, quadric(p, m, (1, 0, 0), const=5))   # x1(x1-1) shifted
    assert len(X2) == 245 and max_on_plane(X2, cube_points(mu2, p)) == 49
