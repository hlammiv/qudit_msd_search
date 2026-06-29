"""Regression tests for the Phase 0 substrate + the M0/M1/K1 findings.

Run:  python -m pytest primepower_msd/tests/test_phase0.py -q   (from the prime_msd repo root)
"""

import numpy as np
import pytest

from primepower_msd.clifford_ring import (
    phase_S,
    fourier_H,
    multiplier_M,
    is_pauli,
    is_clifford,
    is_level3,
    level_of,
)
from primepower_msd.weyl import clock_X, clock_Z, pauli
from primepower_msd.ring import squaring_is_additive, cubing_is_additive
from primepower_msd.single_qudit_gate import monomial_phase_gate, certify_magic
from primepower_msd.hierarchy_search import search, diagonal_clifford, antidifference

DIMS = [2, 4, 8, 16]


# ----------------------------------------------------------------- M0: oracle ground truth (qubit)
def test_qubit_ground_truth():
    """The level oracle must reproduce textbook qubit facts."""
    assert level_of(phase_S(2), 2) == 2                       # S is Clifford
    assert level_of(monomial_phase_gate(2, 2, 1, 8), 2) == 3  # T = diag(1, e^{i pi/4}) is level 3
    assert level_of(fourier_H(2), 2) == 2                     # Hadamard is Clifford


@pytest.mark.parametrize("d", DIMS)
def test_paulis_are_level1(d):
    assert level_of(clock_X(d), d) == 1
    assert level_of(clock_Z(d), d) == 1
    assert is_pauli(pauli(2, 3, d), d)


@pytest.mark.parametrize("d", DIMS)
def test_clifford_generators(d):
    assert is_clifford(phase_S(d), d)        # diag(zeta_{2d}^{x^2}) is the Clifford S, every d
    assert is_clifford(fourier_H(d), d)
    assert is_clifford(multiplier_M(3, d), d)  # 3 is a unit mod 2^k


# ----------------------------------------------------------------- M0: verified algebraic facts
@pytest.mark.parametrize("d", [4, 8, 16, 32])
def test_squaring_cubing_nonadditive(d):
    assert not squaring_is_additive(d)
    assert not cubing_is_additive(d)


@pytest.mark.parametrize("d", [4, 8, 16, 32])
def test_quadratic_gauss_sum(d):
    G = sum(np.exp(2j * np.pi * (x * x) / d) for x in range(d))
    assert np.isclose(G, (1 + 1j) * np.sqrt(d), atol=1e-9)


# ----------------------------------------------------------------- M1: naive monomials are NOT level 3
@pytest.mark.parametrize("d", [4, 8, 16])
def test_naive_monomial_quadratic_not_level3(d):
    """Corrects the literature-map claim: diag(zeta_{4d}^{x^2}) is NOT level 3 for d = 2^k, k>=2."""
    U = monomial_phase_gate(d, 2, 1, 4 * d)
    assert level_of(U, d) >= 4
    cert = certify_magic(U, d)
    assert not cert.is_magic


def test_8th_root_quadratic_is_clifford_on_Z4():
    """diag(zeta_8^{x^2}) on Z_4 is the Clifford S — NOT magic (the corrected precision fact)."""
    U = monomial_phase_gate(4, 2, 1, 8)   # 8 = 2d
    assert is_clifford(U, 4)
    assert np.allclose(U, phase_S(4))


# ----------------------------------------------------------------- K1: existence of single-qudit C_3 magic
@pytest.mark.parametrize("d", DIMS)
def test_k1_existence_count_is_d_times_dminus1(d):
    """Exactly d(d-1) distinct strict-level-3 single-qudit diagonal gates from antidiff(S^a Z^b)."""
    hits = search(d)
    assert len(hits) == d * (d - 1)


@pytest.mark.parametrize("d", [4, 8, 16])
def test_k1_gates_are_genuine_magic_for_powerpower(d):
    """For d = 2^k (k>=2) the level-3 gates pass the anti-collapse certificate (genuine single-ring-qudit)."""
    hits = search(d)
    a, b, sig, U = hits[0]
    cert = certify_magic(U, d)
    assert cert.is_strict_level3 and cert.anticollapse and cert.is_magic
    assert is_level3(U, d) and not is_clifford(U, d)


def test_qubit_level3_is_not_anticollapse_magic():
    """d=2 level-3 gates exist but are NOT flagged 'magic' (shift order 2 -> reducible baseline)."""
    hits = search(2)
    a, b, sig, U = hits[0]
    cert = certify_magic(U, 2)
    assert cert.is_strict_level3
    assert not cert.anticollapse and not cert.is_magic


def test_d4_target_gate_is_wraparound_corrected_quadratic():
    """The d=4 single-ququart magic gate is the wraparound-corrected quadratic exps (0,1,4,13),
    NOT the naive diag(zeta_16^{x^2}) exps (0,1,4,9)."""
    hits = search(4)
    sigs = {h[2] for h in hits}
    assert (0, 1, 4, 13) in sigs          # corrected quadratic IS level-3 magic
    naive = monomial_phase_gate(4, 2, 1, 16)  # exps (0,1,4,9)
    assert level_of(naive, 4) >= 4         # naive is NOT
