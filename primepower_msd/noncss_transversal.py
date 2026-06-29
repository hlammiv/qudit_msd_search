"""V7 machinery (VALIDATED): does a NON-CSS Z_4 stabilizer code admit a TRANSVERSAL single-qudit level-3 gate?

Even-d fix (validated): the raw X^a Z^b has (X^a Z^b)^4 = (-1)^{a.b} I, so when a.b is ODD the matrix
stabilizer group contains -I and there is NO clean codespace.  Use the SIGNED stabilizer
  S'(a,b) = exp(-i*pi*(a.b)/4) X^a Z^b   =>   S'^4 = +I,
and the +1-eigenspace projector (1/4) sum_j S'^j.  On the verified non-CSS [[5,1,3]]_4 code this gives
codespace dim = d^k = 4, idempotent, stabilizers fix it, identity transversal leakage 0.

For a non-CSS code the codewords are not sparse, so we work with the dense codespace projector in
C^(d^n) (feasible for n<=6).  A single-qudit gate G is transversal iff G^{⊗n} preserves the codespace:
(I - Pi) G^{⊗n} Pi == 0.  The induced logical gate is Pi G^{⊗n} Pi in a logical basis; we test whether
it is non-Clifford on the logical qudits.  Validated against the sparse CSS check on a small CSS code.
"""

from __future__ import annotations

import numpy as np

from .weyl import clock_X, clock_Z, omega


def pauli_tensor(a, b, d: int) -> np.ndarray:
    """The Pauli operator (X^{a_i} Z^{b_i})_i as a dense d^n x d^n matrix."""
    n = len(a)
    X, Z = clock_X(d), clock_Z(d)
    Xp = [np.linalg.matrix_power(X, int(ai) % d) for ai in a]
    Zp = [np.linalg.matrix_power(Z, int(bi) % d) for bi in b]
    op = np.array([[1.0 + 0j]])
    for i in range(n):
        op = np.kron(op, Xp[i] @ Zp[i])
    return op


def signed_stabilizer(a, b, d: int) -> np.ndarray:
    """S'(a,b) = exp(-i pi (a.b)/4) X^a Z^b, the order-4 (S'^4=+I) phase-corrected stabilizer (d=4)."""
    ab = int(np.dot(np.asarray(a) % d, np.asarray(b) % d))
    return np.exp(-1j * np.pi * ab / 4.0) * pauli_tensor(a, b, d)


def _order4_proj(S: np.ndarray) -> np.ndarray:
    """Projector onto the +1 eigenspace of an order-4 unitary S: (1/4) sum_{j} S^j."""
    dim = S.shape[0]
    P = np.zeros((dim, dim), dtype=complex)
    M = np.eye(dim, dtype=complex)
    for _ in range(4):
        P += M
        M = M @ S
    return P / 4.0


def codespace_projector(stab_gens, d: int, n: int) -> np.ndarray:
    """Dense projector onto the codespace (+1 eigenspace of the SIGNED stabilizer generators)."""
    dim = d ** n
    Pi = np.eye(dim, dtype=complex)
    for g in stab_gens:
        Pi = Pi @ _order4_proj(signed_stabilizer(g[:n], g[n:], d))
    return Pi


def transversal_test(stab_gens, G: np.ndarray, d: int, n: int):
    """Return (leakage, logical_matrix, dim_codespace) for transversal G^{⊗n} on the code.

    leakage = ||(I-Pi) G^{⊗n} Pi|| (0 => codespace preserved).  logical_matrix = Pi G^{⊗n} Pi restricted
    to an orthonormal codespace basis (None if codespace empty).
    """
    Pi = codespace_projector(stab_gens, d, n)
    # orthonormal codespace basis from columns of Pi
    w, V = np.linalg.eigh((Pi + Pi.conj().T) / 2)
    cols = V[:, w > 0.5]
    kdim = cols.shape[1]
    if kdim == 0:
        return None, None, 0
    Gn = np.array([[1.0 + 0j]])
    for _ in range(n):
        Gn = np.kron(Gn, G)
    GPi_cols = Gn @ cols
    # leakage = component outside codespace
    leak = np.linalg.norm(GPi_cols - cols @ (cols.conj().T @ GPi_cols))
    logical = cols.conj().T @ GPi_cols      # kdim x kdim logical matrix
    return float(leak), logical, kdim
