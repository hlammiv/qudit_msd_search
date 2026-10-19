"""Universality (sector-escape) criterion — the table in UNIVERSALITY_CHECK.md, symbolic == dense."""

import itertools as it

import numpy as np
import pytest

from primepower_msd import a2
from primepower_msd import phase_tables as PT


def test_known_gates_table():
    preserve = {
        "S (Clifford)": PT.clifford_table(1, 0),
        "Z": PT.pauli_table([1]),
        "CZ i^{xy}": PT.cz_layer(2, [(0, 1, 1)]),
        "ring-CCZ i^{xyz}": PT.ccz_layer(3, [(0, 1, 2, 1)]),
        "i^{x^3}": np.array([0, 4, 0, 12]),
    }
    escape = {
        "carrier zeta16^[0,1,4,13]": np.array(PT.CARRIER_Q),
        "zeta8^{x^3}": np.array([0, 2, 0, 6]),
        "zeta16^{x}": np.array([0, 1, 2, 3]),
    }
    for name, P in preserve.items():
        assert not a2.escapes_A2(P), name
        assert not a2.escapes_A2_dense(PT.to_unitary(P), P.ndim), name
    for name, P in escape.items():
        assert a2.escapes_A2(P), name
        assert a2.escapes_A2_dense(PT.to_unitary(P), P.ndim), name


def test_mod4_polynomial_phases_preserve_A2():
    rng = np.random.default_rng(7)
    for m in (1, 2, 3):
        mons = [c for deg in (1, 2, 3) for c in it.combinations_with_replacement(range(m), deg)]
        for _ in range(40):
            coeffs = rng.integers(0, 4, size=len(mons))
            P = PT.from_poly_mod4(coeffs, mons, m)
            assert not a2.escapes_A2(P)
            if m <= 2:
                assert not a2.escapes_A2_dense(PT.to_unitary(P), m)


def test_symbolic_equals_dense_single_qudit_all_tables():
    bad = 0
    for rest in it.product(range(16), repeat=3):
        P = np.array((0,) + rest)
        if a2.escapes_A2(P) != a2.escapes_A2_dense(PT.to_unitary(P), 1):
            bad += 1
    assert bad == 0


def test_symbolic_equals_dense_two_qudits_random():
    rng = np.random.default_rng(11)
    for _ in range(120):
        P = PT.normalize(rng.integers(0, 16, size=(4, 4)))
        assert a2.escapes_A2(P) == a2.escapes_A2_dense(PT.to_unitary(P), 2)


def test_d8_sectors_symbolic_equals_dense():
    rng = np.random.default_rng(5)
    for _ in range(150):
        P = PT.normalize(rng.integers(0, 32, size=8), 8)
        U = PT.to_unitary(P, 8)
        for j in (1, 2):
            assert a2.preserves_sector(P, j, 8) == a2.preserves_sector_dense(U, 1, 8, j)


def test_escape_witness_consistent():
    assert a2.escape_witness(np.array(PT.CARRIER_Q)) is not None
    assert a2.escape_witness(PT.clifford_table(3, 1)) is None


def test_entangling_escaper_exists_and_is_level3():
    """The one-body 'lemma' is FALSE on tables: this 2-ququart gate escapes A2 through its
    entangling part (zero one-body marginals), is level 3 and has order 8."""
    P = np.array([[0, 0, 0, 0], [0, 2, 0, 10], [0, 0, 8, 8], [0, 10, 8, 10]])
    assert PT.level(P) == 3
    assert PT.level_of_dense(PT.to_unitary(P), 2) == 3
    assert a2.escapes_A2(P) and a2.escapes_A2_dense(PT.to_unitary(P), 2)
    assert all(not np.any(q) for q in PT.qudit_marginals(P))
    assert a2.escapes_A2(PT.mixed_part(P))
