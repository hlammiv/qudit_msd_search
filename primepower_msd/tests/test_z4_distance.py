"""Z_4 meet-in-the-middle distance: exact against brute force on stored and random codes."""

import random
from pathlib import Path

import numpy as np

from primepower_msd import z4_distance as ZD
from primepower_msd.ringlinalg import right_kernel_howell
from primepower_msd.symplectic import Z4Code, brute_distance, symp

PKG = Path(__file__).resolve().parents[1]
CODE513 = [[1, 2, 1, 1, 2, 1, 0, 3, 3, 1], [1, 1, 1, 2, 1, 2, 1, 3, 0, 2],
           [3, 2, 0, 1, 2, 0, 1, 1, 1, 0], [2, 1, 2, 0, 2, 2, 1, 3, 1, 3]]


def test_stored_codes():
    for name, exp in [("witness_8_2_2_hosts_carrier.npy", 2), ("code_10_4_3_gamma0834.npy", 3),
                      ("code_12_5_3_gamma0797.npy", 3), ("code_12_2_3_d3_k2.npy", 3)]:
        c = Z4Code.load(PKG / name)
        w, cert = ZD.symplectic_distance(c.rows, c.n, d_max=4)
        assert (w, cert) == (exp, True)
        dmin, cert2, _ = ZD.code_distance(c, d_max=4)
        assert dmin == exp and cert2
    assert ZD.symplectic_distance(CODE513, 5, d_max=4) == (3, True)


def test_random_noncss_matches_brute():
    rng = random.Random(2)
    for _ in range(40):
        n = rng.choice([4, 5, 6])
        rows = []
        target = rng.randrange(1, n)
        while len(rows) < target:
            v = [rng.randrange(4) for _ in range(2 * n)]
            if any(v) and all(symp(v, g, n) == 0 for g in rows):
                rows.append(v)
        bd = brute_distance(rows, n, cap=3)
        w, _ = ZD.symplectic_distance(rows, n, d_max=3)
        assert bd == (w if w else 4)


def test_random_css_matches_brute_and_handles_stabilizer_collisions():
    rng = random.Random(3)
    checked = 0
    for _ in range(60):
        n = rng.choice([4, 5, 6])
        MX = [[rng.randrange(4) for _ in range(n)] for _ in range(rng.randrange(1, 3))]
        K = right_kernel_howell(np.array(MX), 4)
        MZ = [list(map(int, K[i])) for i in range(len(K)) if rng.random() < 0.5]
        try:
            c = Z4Code.from_css(MX, MZ, n)
        except AssertionError:
            continue
        dmin, dX, dZ, cert = ZD.css_distance(MX, MZ, n, d_max=3)
        bd = brute_distance(c.rows, n, cap=3)
        assert bd == (dmin if dmin else 4)
        checked += 1
    assert checked >= 30


def test_lightest_dual_element_is_a_stabilizer_no_false_exit():
    """MX = [2,2,0,0]: the weight-1 vector (2,0,0,0) is in MX^perp (2*2=0) and weight-1 Z (0,0,1,0) too;
    with MZ containing (2,0,0,0) and (0,0,1,0),(0,0,0,1) the lightest Z-logical must skip them."""
    n = 4
    MX = [[2, 2, 0, 0]]
    MZ = [[2, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
    c = Z4Code.from_css(MX, MZ, n)
    dmin, dX, dZ, cert = ZD.css_distance(MX, MZ, n, d_max=3)
    assert dZ == brute_distance(c.rows, n, cap=3) or dmin == brute_distance(c.rows, n, cap=3)
    assert dmin == brute_distance(c.rows, n, cap=3)


def test_memory_guard():
    import pytest
    with pytest.raises(MemoryError):
        ZD.weight_w_logical(np.zeros((1, 60), dtype=int), [], 60, 6, "noncss", 4, max_left_entries=1000)
