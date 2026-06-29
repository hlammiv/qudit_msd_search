"""Cyclic-clock Heisenberg-Weyl (generalized Pauli) group for a single d-level qudit.

The carrier is the *cyclic* clock group  X|j> = |j+1 mod d>,  Z|j> = omega_d^j |j>  with
ZX = omega_d XZ.  For the prime-power target this is the MANDATORY substrate: the standard
tensor/field Pauli group on p^r splits into r prime-p qudits (Cui-Gottesman-Krishna,
arXiv:1608.06596) — the very collapse the project must avoid — whereas the cyclic clock group
is indecomposable for d = 2^k (the shift X has order d, an order no tensor product of smaller
2-power qudits can match).  Everything is an exact d x d complex matrix.
"""

from __future__ import annotations

import numpy as np


def omega(d: int) -> complex:
    """Primitive d-th root of unity exp(2*pi*i/d)."""
    return np.exp(2j * np.pi / d)


def clock_X(d: int) -> np.ndarray:
    """Shift operator X|j> = |j+1 mod d> (column j has its 1 in row (j+1) mod d)."""
    X = np.zeros((d, d), dtype=complex)
    for j in range(d):
        X[(j + 1) % d, j] = 1.0
    return X


def clock_Z(d: int) -> np.ndarray:
    """Clock operator Z|j> = omega_d^j |j>."""
    w = omega(d)
    return np.diag([w ** j for j in range(d)])


def Xpow(a: int, d: int) -> np.ndarray:
    """X^a as an explicit shift-by-a matrix (a taken mod d)."""
    a %= d
    X = np.zeros((d, d), dtype=complex)
    for j in range(d):
        X[(j + a) % d, j] = 1.0
    return X


def Zpow(b: int, d: int) -> np.ndarray:
    """Z^b = diag(omega_d^{b j})."""
    w = omega(d)
    return np.diag([w ** ((b * j) % d) for j in range(d)])


def pauli(a: int, b: int, d: int) -> np.ndarray:
    """Generalized Pauli X^a Z^b:  (X^a Z^b)|j> = omega_d^{b j} |j + a mod d>."""
    w = omega(d)
    P = np.zeros((d, d), dtype=complex)
    for j in range(d):
        P[(j + a) % d, j] = w ** ((b * j) % d)
    return P


def pauli_shift_order(d: int) -> int:
    """Order of the shift generator X in the clock group (= d).

    Anti-collapse certificate (substrate level): a Clifford map is a Pauli-group automorphism
    and preserves element orders.  For d = 2^k (k >= 2) the cyclic clock group contains the
    order-d shift X, while any tensor product of smaller 2-power qudits has shift-type generators
    of order at most the largest factor < d.  Hence the cyclic Z_{2^k} Pauli group is NOT
    isomorphic to a tensor product of smaller-qudit Pauli groups — the resource cannot be
    relabelled into "qubit MSD in GF(4) disguise" at the Pauli-group level.
    """
    return d
