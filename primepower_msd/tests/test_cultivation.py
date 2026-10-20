import itertools as it
from math import comb

import numpy as np
import pytest

from primepower_msd import cultivation as CU
from primepower_msd import phase_tables as PT
from primepower_msd.a2 import escapes_all_sectors

pytest.importorskip("stim")


@pytest.mark.parametrize("d", [4, 8])
def test_escape_iff_odd_check_exponent(d):
    M = 4 * d
    seen = set()
    for A, B, C in it.product(range(M), repeat=3):
        q = np.array([(A * y + B * comb(y, 2) + C * comb(y, 3)) % M for y in range(2 * d)])
        if not np.array_equal(q[:d], q[d:]) or tuple(q[:d]) in seen:
            continue
        q = q[:d]
        seen.add(tuple(q))
        if PT.level(q, d) != 3:
            continue
        g = CU.check_table(q, d)
        assert PT.level(g, d) <= 2
        assert (CU.s_exponent(g, d) % 2 == 1) == escapes_all_sectors(q, d)


@pytest.mark.parametrize("d,dist", [(4, 3), (4, 5), (8, 3), (8, 5), (16, 3)])
def test_transversal_t_check_on_color_code(d, dist):
    r = CU.verify_check(dist, d)
    assert r["n"] == (3 * dist * dist + 1) // 4 and r["k"] == 1
    assert r["hosted"] and r["x_all_ones_logical"]
    assert r["logical_s_exponent"] % 2 == 1
    diff = {(a - b) % (4 * d) for a, b in zip(r["logical_g"], r["physical_g"])}
    assert len(diff) == 1                              # logical check = T_L X_L T_L^dag (same carrier, up to phase)


def test_dense_steane_z4():
    d, M = 4, 16
    code, tau = CU.signed_color_code(3, d)
    n, q = code.n, CU.escaping_carrier(d)
    z = np.exp(2j * np.pi / M)
    span = {tuple([0] * n)}
    for g0 in np.array(code.x_gens) % d:
        span = {tuple((np.array(s) + k * g0) % d) for s in span for k in range(d)}
    idx = lambda x: int(sum(int(v) * d ** i for i, v in enumerate(x)))
    L = []
    for j in range(d):
        v = np.zeros(d ** n, dtype=complex)
        for s in span:
            v[idx((np.array(s) + j) % d)] += 1
        L.append(v / np.linalg.norm(v))
    L = np.array(L)

    def apply_check(v):
        out = np.zeros_like(v)
        for k in np.nonzero(np.abs(v) > 1e-12)[0]:
            x = np.array([(k // d ** i) % d for i in range(n)])
            out[idx((x + 1) % d)] += v[k] * z ** sum(int(t) * int(q[(xi + 1) % d] - q[xi]) for t, xi in zip(tau, x))
        return out

    CL = np.array([[L[a].conj() @ apply_check(L[b]) for b in range(d)] for a in range(d)])
    assert max(np.linalg.norm(apply_check(L[b]) - L.T @ CL[:, b]) for b in range(d)) < 1e-9
    # |T_L> = T_L |+_L> is a +1 eigenstate of the transversal check
    TL = np.diag(z ** q)
    tl = TL @ np.ones(d) / 2
    assert np.allclose(CL @ tl, tl)
