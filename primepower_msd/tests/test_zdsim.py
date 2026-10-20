import numpy as np
import pytest

from primepower_msd.zdsim import Circuit, FrameSim, check_determinism


def pauli(x, z, d):
    X = np.roll(np.eye(d), 1, axis=0)
    Z = np.diag(np.exp(2j * np.pi * np.arange(d) / d))
    return np.linalg.matrix_power(X, x % d) @ np.linalg.matrix_power(Z, z % d)


def same_up_to_phase(A, B):
    k = np.argmax(np.abs(B))
    ph = A.flat[k] / B.flat[k]
    return abs(abs(ph) - 1) < 1e-9 and np.allclose(A, ph * B)


@pytest.mark.parametrize("d", [2, 3, 4, 8])
def test_single_qudit_rules(d):
    w = np.exp(2j * np.pi / d)
    Fm = np.array([[w ** (j * k) for j in range(d)] for k in range(d)]) / np.sqrt(d)
    Sm = np.diag(np.exp(2j * np.pi * np.arange(d) ** 2 / (2 * d)))
    for a in range(1, d):
        Sa = np.linalg.matrix_power(Sm, a)
        for x in range(d):
            for z in range(d):
                P = pauli(x, z, d)
                assert same_up_to_phase(Fm @ P @ Fm.conj().T, pauli(-z, x, d))
                assert same_up_to_phase(Fm.conj().T @ P @ Fm, pauli(z, -x, d))
                if d % 2 == 0:                                   # S = diag(zeta_{2d}^{x^2}) needs even d
                    assert same_up_to_phase(Sa @ P @ Sa.conj().T, pauli(x, z + a * x, d))


@pytest.mark.parametrize("d", [2, 4])
def test_csum_rule(d):
    U = np.zeros((d * d, d * d))
    for a in range(d):
        for b in range(d):
            U[a * d + (b + 2 * a) % d, a * d + b] = 1          # w = 2: |a,b> -> |a, b+2a>
    for xc in range(d):
        for zc in range(d):
            for xt in range(d):
                for zt in range(d):
                    P = np.kron(pauli(xc, zc, d), pauli(xt, zt, d))
                    Q = np.kron(pauli(xc, zc - 2 * zt, d), pauli(xt + 2 * xc, zt, d))
                    assert same_up_to_phase(U @ P @ U.T, Q)


def test_frame_matches_rules_and_determinism():
    c = Circuit(4); a, b = c.q(2)
    c.R([a, b]); c.F([a]); c.CSUM([(a, b)]); c.CSUM([(a, b)], w=-1); c.FD([a])
    m = c.M([a, b]); c.detector([(m[0], 1)]); c.detector([(m[1], 1)])
    assert check_determinism(c) == []
    c2 = Circuit(4); a, b = c2.q(2)
    c2.R([a, b]); c2.F([a]); c2.CSUM([(a, b)]); m = c2.M([a, b]); c2.detector([(m[1], 1)])
    assert check_determinism(c2) != []            # a single GHZ half is random
    c2.detector([(m[0], 1), (m[1], -1)])
    assert len(check_determinism(c2)) == 1         # the parity is deterministic
