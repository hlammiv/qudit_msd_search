"""V7 VALIDATED dense codespace + transversal-gate tester for NON-CSS Z_4 (even-d) stabilizer codes.

This module FIXES the even-d Weyl phase bug that made the old `noncss_transversal.codespace_projector`
return the wrong rank.

THE BUG.  A stabilizer generator on n qudits is g = prod_i (X^{a_i} Z^{b_i}).  Using the naive Pauli
X^a Z^b, for EVEN d each single-qudit factor obeys (X^a Z^b)^d = (-1)^{ab} I, hence the tensor obeys
    g^d = (-1)^{sum_i a_i b_i} I .
When sum_i a_i b_i is ODD, g^d = -I, so g has 2d-th-root eigenvalues, +1 is NOT an eigenvalue, and the
projector (1/d) sum_{j=0}^{d-1} g^j collapses to (near) rank 0 on that factor -> wrong codespace dim
(the observed "rank 4 instead of 16").

THE FIX (displacement / tau convention).  Use D(a,b) = tau^{ab} X^a Z^b with tau = e^{i*pi*(d+1)/d}
(a primitive 2d-th root: tau^2 = omega_d, tau^d = -1 for even d).  Then a single-qudit displacement obeys
D(a,b)^d = I exactly, so for the tensor generator g = prod_i D(a_i,b_i) we get g^d = I ALWAYS, and
P_g = (1/d) sum_{j=0}^{d-1} g^j is a true orthogonal projector onto g's +1 eigenspace.  For commuting
generators the P_g commute and Pi = prod_g P_g is the joint codespace projector of rank d^{n-k}.

API:
    displacement(a, b, d)            single-qudit D(a,b)
    pauli_tensor(a_vec, b_vec, d)    n-qudit displacement (tau convention)
    codespace_projector(rows, d, n)  dense joint +1 projector (rows = list of (a|b) in Z_d^{2n})
    codespace_basis(Pi)              orthonormal columns spanning the codespace
    transversal_test(rows, G, d, n)  -> (leakage, logical_matrix, kdim)

Feasible for n <= 6 at d = 4 (d^6 = 4096 -> 4096x4096 complex ~ 256 MB ceiling).
"""

from __future__ import annotations

import numpy as np

from .weyl import clock_X, clock_Z, omega

_TOL = 1e-9


def tau(d: int) -> complex:
    """Primitive 2d-th root tau = e^{i*pi*(d+1)/d};  tau^2 = omega_d,  tau^d = (-1)^d (= -1 for even d)."""
    return np.exp(1j * np.pi * (d + 1) / d)


def displacement(a: int, b: int, d: int) -> np.ndarray:
    """Single-qudit displacement D(a,b) = tau^{ab} X^a Z^b;  satisfies D(a,b)^d = I for all d."""
    a %= d
    b %= d
    X = clock_X(d)
    Z = clock_Z(d)
    Xa = np.linalg.matrix_power(X, a)
    Zb = np.linalg.matrix_power(Z, b)
    return (tau(d) ** ((a * b) % (2 * d))) * (Xa @ Zb)


def pauli_tensor(a_vec, b_vec, d: int) -> np.ndarray:
    """n-qudit displacement g = prod_i D(a_i, b_i) as a dense d^n x d^n matrix (tau convention)."""
    n = len(a_vec)
    op = np.array([[1.0 + 0j]])
    for i in range(n):
        op = np.kron(op, displacement(int(a_vec[i]), int(b_vec[i]), d))
    return op


def single_proj(g: np.ndarray, d: int) -> np.ndarray:
    """Projector (1/d) sum_{j=0}^{d-1} g^j onto the +1 eigenspace of g (valid because g^d = I)."""
    dim = g.shape[0]
    P = np.zeros((dim, dim), dtype=complex)
    M = np.eye(dim, dtype=complex)
    for _ in range(d):
        P += M
        M = M @ g
    return P / d


def codespace_projector(rows, d: int, n: int) -> np.ndarray:
    """Dense joint +1 projector for stabilizer generators given as rows (a|b) of length 2n over Z_d.

    Assumes the rows pairwise commute (symplectic form a.b' - a'.b == 0 mod d); the caller should
    have checked this (the centralizer/code-validity is a precondition of being a stabilizer code).
    """
    Pi = np.eye(d ** n, dtype=complex)
    for row in rows:
        a, b = row[:n], row[n:]
        g = pauli_tensor(a, b, d)
        Pi = Pi @ single_proj(g, d)
    return Pi


def codespace_basis(Pi: np.ndarray, tol: float = 1e-7) -> np.ndarray:
    """Orthonormal basis (columns) for the range of the Hermitian projector Pi via its top eigenvectors."""
    H = (Pi + Pi.conj().T) / 2.0
    w, V = np.linalg.eigh(H)
    cols = V[:, w > 0.5]
    return cols


def transversal_gate(G: np.ndarray, n: int) -> np.ndarray:
    """G^{⊗n} as a dense d^n x d^n matrix."""
    Gn = np.array([[1.0 + 0j]])
    for _ in range(n):
        Gn = np.kron(Gn, G)
    return Gn


def transversal_test(rows, G: np.ndarray, d: int, n: int, Pi: np.ndarray | None = None):
    """(leakage, logical_matrix, kdim) for transversal G^{⊗n} on the code defined by `rows`.

    leakage = || (I - Pi) G^{⊗n} Pi ||_F  (0  <=>  codespace preserved).
    logical_matrix = B^dag G^{⊗n} B in an orthonormal codespace basis B (kdim x kdim), or None if empty.
    """
    if Pi is None:
        Pi = codespace_projector(rows, d, n)
    B = codespace_basis(Pi)
    kdim = B.shape[1]
    if kdim == 0:
        return None, None, 0
    Gn = transversal_gate(G, n)
    GB = Gn @ B
    leak = float(np.linalg.norm(GB - B @ (B.conj().T @ GB)))
    logical = B.conj().T @ GB
    return leak, logical, kdim
