"""Single-qudit Clifford group over the cyclic clock group, and an exact Clifford-hierarchy
level oracle.

Membership tests are *operational* and exact (matrix computation), so the even-d phase-convention
subtleties (the 2d-doubling) are handled automatically rather than by hand:

  * Pauli group  C_1 = { phase . X^a Z^b }.
  * is_clifford(U)  <=>  U X U^dag and U Z U^dag are both Paulis (X, Z generate C_1, and
    conjugation is a homomorphism).
  * is_level3(U)    <=>  U X U^dag and U Z U^dag are both Clifford  (i.e. U in C_3).
  * level_of(U)     in {1, 2, 3, >=4}.

This is the de Silva-Lautsch / Cui-Gottesman-Krishna single-qudit hierarchy made concrete for the
cyclic ring Z_{2^k}, which those works leave open beyond prime d.
"""

from __future__ import annotations

import numpy as np

from .weyl import omega, clock_X, clock_Z, Xpow, Zpow, pauli

_TOL = 1e-9


# --------------------------------------------------------------------------- Clifford generators
def phase_S(d: int) -> np.ndarray:
    """Canonical Clifford phase gate S = diag(mu^{j^2}), mu = exp(2*pi*i/(2d)) (even-d doubling).

    S X S^dag = mu . Z X (a Pauli), so S is Clifford.  On Z_4, mu = zeta_8, i.e. the 8th-root
    quadratic diag(zeta_8^{x^2}) is exactly this Clifford S — NOT a magic gate (the corrected fact
    from the literature survey).
    """
    mu = np.exp(2j * np.pi / (2 * d))
    return np.diag([mu ** ((j * j) % (2 * d)) for j in range(d)])


def fourier_H(d: int) -> np.ndarray:
    """Fourier (Hadamard) gate H[j,k] = omega_d^{jk}/sqrt(d); Clifford for all d."""
    w = omega(d)
    H = np.array([[w ** ((j * k) % d) for k in range(d)] for j in range(d)], dtype=complex)
    return H / np.sqrt(d)


def multiplier_M(a: int, d: int) -> np.ndarray:
    """Multiplier M_a|j> = |a j mod d>, a a unit of Z_d; Clifford."""
    M = np.zeros((d, d), dtype=complex)
    for j in range(d):
        M[(a * j) % d, j] = 1.0
    return M


# --------------------------------------------------------------------------- equality / structure
def _equal_up_to_phase(A: np.ndarray, B: np.ndarray, tol: float = _TOL) -> bool:
    """True iff A = e^{i phi} B, i.e. A^dag B is a unit-modulus scalar multiple of the identity."""
    M = A.conj().T @ B
    diag = np.diag(M)
    if not np.allclose(M - np.diag(diag), 0.0, atol=tol):
        return False
    if not np.allclose(diag, diag[0], atol=tol):
        return False
    return abs(abs(diag[0]) - 1.0) < tol


def _is_monomial(V: np.ndarray, tol: float = _TOL) -> bool:
    """True iff V has exactly one nonzero entry per row and per column."""
    mask = np.abs(V) > tol
    return bool(np.all(mask.sum(axis=0) == 1) and np.all(mask.sum(axis=1) == 1))


# --------------------------------------------------------------------------- hierarchy membership
def is_pauli(V: np.ndarray, d: int, tol: float = _TOL) -> bool:
    """True iff V = (global phase) . X^a Z^b for some a, b in Z_d."""
    if not _is_monomial(V, tol):
        return False
    # shift a: in column 0 the nonzero sits in row a (since (X^a Z^b)|0> ∝ |a>).
    col0 = np.where(np.abs(V[:, 0]) > tol)[0]
    if len(col0) != 1:
        return False
    a = int(col0[0])
    D = V @ Xpow(a, d).conj().T  # remove the shift; should be diagonal = phase . Z^b
    if not np.allclose(D - np.diag(np.diag(D)), 0.0, atol=tol):
        return False
    dvals = np.diag(D)
    c = dvals[0]
    if abs(abs(c) - 1.0) > tol and abs(c) < tol:
        return False
    ratios = dvals / c
    w = omega(d)
    for b in range(d):
        if np.allclose(ratios, [w ** ((b * j) % d) for j in range(d)], atol=1e-7):
            return abs(abs(c) - 1.0) < 1e-7
    return False


def is_clifford(U: np.ndarray, d: int, tol: float = _TOL) -> bool:
    """True iff U is in the Clifford group C_2 (normalizes the cyclic Pauli group)."""
    Ud = U.conj().T
    X, Z = clock_X(d), clock_Z(d)
    return is_pauli(U @ X @ Ud, d, tol) and is_pauli(U @ Z @ Ud, d, tol)


def is_level3(U: np.ndarray, d: int, tol: float = _TOL) -> bool:
    """True iff U is in C_3 (conjugates every Pauli into the Clifford group)."""
    Ud = U.conj().T
    X, Z = clock_X(d), clock_Z(d)
    return is_clifford(U @ X @ Ud, d, tol) and is_clifford(U @ Z @ Ud, d, tol)


def level_of(U: np.ndarray, d: int, tol: float = _TOL) -> int:
    """Clifford-hierarchy level of U: 1 (Pauli), 2 (Clifford), 3 (C_3), or 4 meaning '>= 4 / not in C_3'."""
    if is_pauli(U, d, tol):
        return 1
    if is_clifford(U, d, tol):
        return 2
    if is_level3(U, d, tol):
        return 3
    return 4
