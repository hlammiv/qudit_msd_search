"""Cyclic triorthogonal oracle builder (arXiv:2610.08012)."""
import numpy as np

from qmsd.cyclic import (admissible_A, admissible_B, construction_A, construction_B,
                         cyclotomic_cosets, sumset)
from qmsd.performance import profile_exact
from qmsd.triorthogonal import is_triorthogonal_cubic, is_cube_mod_p, logical_lambdas

# Example 14: explicit H of [[19,1,3]]_3, logical row first.
EX14 = np.array([
    [0, 0, 0, 0, 0, 2, 2, 1, 1, 2, 0, 2, 0, 1, 0, 1, 2, 2, 1],
    [1, 0, 0, 0, 0, 1, 0, 1, 0, 2, 2, 1, 2, 2, 1, 2, 2, 0, 1],
    [0, 1, 0, 0, 0, 1, 2, 2, 0, 1, 2, 0, 1, 1, 2, 0, 0, 0, 2],
    [0, 0, 1, 0, 0, 2, 0, 0, 0, 2, 1, 1, 0, 2, 1, 0, 2, 2, 1],
    [0, 0, 0, 1, 0, 1, 0, 2, 2, 1, 2, 2, 1, 2, 2, 0, 1, 0, 1],
    [0, 0, 0, 0, 1, 1, 2, 2, 1, 0, 1, 0, 2, 0, 2, 1, 1, 2, 2],
])


def test_example3_cosets_and_sumsets():
    assert cyclotomic_cosets(3, 8) == [(0,), (1, 3), (2, 6), (4,), (5, 7)]
    assert sumset({1, 3}, 2, 8) == {2, 4, 6}
    assert sumset({1, 3}, 3, 8) == {1, 3, 5, 7}


def test_example14_not_self_orthogonal_but_B_admissible():
    S = {2, 5, 6, 14, 15, 18}
    assert 0 in sumset(S, 2, 20) and admissible_B(S, 20) and not admissible_A(S, 20)


def test_example14_matches_construction_B():
    assert is_triorthogonal_cubic(EX14, 1, 3)
    G0, H1 = construction_B(3, 20, {2, 5, 6, 14, 15, 18})
    built = profile_exact(G0, H1, 3)
    given = profile_exact(EX14[1:], EX14[:1], 3)
    assert built.accept.D == given.accept.D
    assert built.logical["block"].D == given.logical["block"].D


def test_lambdas_and_cube_classes():
    for p, n, S, cons, lam in [(7, 8, {1, 4, 7}, construction_B, 6),
                               (7, 12, {1, 2, 7}, construction_A, 5),
                               (3, 40, {1, 2, 3, 5, 6, 9, 14, 15, 18, 27}, construction_A, 1)]:
        G0, H1 = cons(p, n, S)
        H = np.vstack([H1, G0])
        assert is_triorthogonal_cubic(H, 1, p)
        assert logical_lambdas(H, 1, p) == [lam]
    assert is_cube_mod_p(6, 7) and not is_cube_mod_p(5, 7)  # [[12,1,3]]_7 is the odd one out
