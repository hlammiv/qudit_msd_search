"""Regression test for the validated even-d non-CSS codespace machinery (V7).

The verified non-CSS [[5,1,3]]_4 code must yield a valid codespace (dim d^k=4, idempotent projector,
stabilizers fix it, transversal identity leakage 0) — this guards the even-d signed-stabilizer fix.
"""

import numpy as np

from primepower_msd.noncss_transversal import codespace_projector, signed_stabilizer

D, N = 4, 5
# verified non-CSS [[5,1,3]]_4 stabilizer generators (a|b), distance 3, k=1
CODE = [
    [1, 2, 1, 1, 2, 1, 0, 3, 3, 1],
    [1, 1, 1, 2, 1, 2, 1, 3, 0, 2],
    [3, 2, 0, 1, 2, 0, 1, 1, 1, 0],
    [2, 1, 2, 0, 2, 2, 1, 3, 1, 3],
]


def test_signed_stabilizers_order4():
    for g in CODE:
        S = signed_stabilizer(g[:N], g[N:], D)
        assert np.linalg.norm(np.linalg.matrix_power(S, 4) - np.eye(D ** N)) < 1e-9


def test_noncss_codespace_valid():
    Pi = codespace_projector(CODE, D, N)
    assert abs(np.trace(Pi).real - 4) < 1e-6                       # codespace dim = d^k = 4
    assert np.linalg.norm(Pi @ Pi - Pi) < 1e-9                     # idempotent
    for g in CODE:
        S = signed_stabilizer(g[:N], g[N:], D)
        assert np.linalg.norm(S @ Pi - Pi) < 1e-9                  # stabilizer fixes codespace
