import itertools as it

import numpy as np
import pytest

from primepower_msd import rgt
from primepower_msd.cultivation import check_layer, check_table, escaping_carrier, s_exponent
from primepower_msd.meas_distill import (analyze, balance_report, gamma_hhpw, logical_x_basis, normal_basis,
                                         signed_rm)
from primepower_msd.phase_tables import level, modulus
from primepower_msd.symplectic import Z4Code


@pytest.mark.parametrize("d", [4, 8, 16])
def test_t_inverse_costs_one_t(d):
    q = escaping_carrier(d); x = np.arange(d)
    assert level((q[(-x) % d] + q) % modulus(d), d) <= 2          # T^{-1} = N T N^dag (Clifford)


@pytest.mark.parametrize("r,m", [(1, 3), (1, 4), (2, 5), (2, 6), (3, 7)])
def test_signed_rm_exact_balance(r, m):
    S, tau = signed_rm(r, m)
    assert balance_report(S, tau)["exact"]


@pytest.mark.parametrize("r,m,d", [(1, 3, 4), (1, 4, 4), (1, 4, 8), (2, 5, 4)])
def test_all_ones_check_hosted_and_separable(r, m, d):
    S, tau = signed_rm(r, m); n = S.shape[1]; M = modulus(d)
    code = Z4Code.from_css(S % d, (S * tau) % d, n, d=d)
    q = escaping_carrier(d)
    assert rgt.hosts(code, check_layer(code, tau, q), "support").hosted
    N = normal_basis(logical_x_basis(S, tau, d), tau, d)
    assert N is not None and len(N) == n - 2 * len(S)
    g = check_table(q, d)
    Phi = lambda v: int(sum(int(t) * g[int(a) % d] for t, a in zip(tau, v))) % M
    k = len(N)
    single = [[Phi((x * N[j]) % d) for x in range(d)] for j in range(k)]
    rng = np.random.default_rng(0)
    for _ in range(300):
        xs = rng.integers(0, d, k)
        v = (xs @ N) % d
        assert (Phi(v) - sum(single[j][xs[j]] for j in range(k)) + (k - 1) * Phi(np.zeros(n))) % M == 0
    for j in range(k):
        assert s_exponent(np.array(single[j]), d) % 2 == 1                      # escaping per logical
        assert int((N[j] * tau).sum()) % 2 == 1                                 # odd X exponent c_j


def test_distance_and_gamma_63():
    S, tau = signed_rm(2, 6)
    a = analyze(S, tau, 4, D_max=6)
    assert a["k"] == 21 and a["normal"] and a["d_Z"] == ">6"                   # D = 7 (binary upper bound)
    assert abs(gamma_hhpw(63, 21, 7) - np.log(19) / np.log(7)) < 1e-12
