"""V1 (FIELD) machinery: GF(4) single-qudit Weyl-Heisenberg substrate, codespace tester,
transversal level-3 diagonal gate, and the ANTI-COLLAPSE certificate.

A "single GF(4)-qudit" is C^4 with basis indexed by GF(4) = {0, 1, w, w^2}.  Its Weyl group
uses the ADDITIVE group of GF(4), which is (Z_2)^2 (characteristic 2): every additive shift
X_a satisfies X_a^2 = X_{2a} = X_0 = I, so EVERY shift has order 2 -- NOT 4.  This is the
field/2-qubit collapse the project rejects: the GF(4) Pauli group is literally the 2-qubit
Pauli group, and the candidate magic gate diag(1,1,1,i) is the entangling 2-qubit
controlled-S, not a genuine single-(cyclic-Z_4)-qudit gate.

This module BUILDS the resource honestly so the collapse can be measured rather than asserted.
"""

from __future__ import annotations

import itertools

import numpy as np
import galois

GF4 = galois.GF(4)

# Canonical ordering of the 4 field elements as basis labels |0>,|1>,|w>,|w^2>.
# galois represents GF(4) elements as ints 0,1,2,3 with 2 == w (primitive), 3 == w+1 == w^2.
ELEMS = [GF4(0), GF4(1), GF4(2), GF4(3)]          # 0, 1, w, w^2
IDX = {int(e): i for i, e in enumerate(ELEMS)}     # int-value -> basis position


# --------------------------------------------------------------------------- field trace
def trace(x) -> int:
    """Field trace Tr: GF(4) -> GF(2), Tr(x) = x + x^2, returned as int 0/1.

    Tr(0)=0, Tr(1)=0, Tr(w)=1, Tr(w^2)=1.  The additive character of the GF(4) Weyl group is
    (-1)^{Tr(.)}; characteristic 2 is exactly why the shift order collapses to 2.
    """
    x = GF4(x)
    return int(x + x ** 2)  # element of GF(2) embedded as 0/1


# --------------------------------------------------------------------------- single-qudit Weyl
def weyl_X(a) -> np.ndarray:
    """Additive shift X_a|x> = |x + a>, a in GF(4).  Real permutation matrix, X_a^2 = I."""
    a = GF4(a)
    M = np.zeros((4, 4), dtype=complex)
    for i, x in enumerate(ELEMS):
        j = IDX[int(x + a)]
        M[j, i] = 1.0
    return M


def weyl_Z(b) -> np.ndarray:
    """Phase clock Z_b|x> = (-1)^{Tr(b x)} |x>, b in GF(4).  Real diagonal +-1, Z_b^2 = I."""
    b = GF4(b)
    return np.diag([(-1.0) ** trace(b * x) for x in ELEMS]).astype(complex)


def weyl(a, b) -> np.ndarray:
    """W(a,b) = X_a Z_b (single GF(4) qudit)."""
    return weyl_X(a) @ weyl_Z(b)


# --------------------------------------------------------------------------- commutation / form
def trace_symplectic(a, b, ap, bp) -> int:
    """Trace-symplectic form <(a|b),(a'|b')> = Tr(sum_i (a_i b'_i + a'_i b_i)) in GF(2) (char 2).

    Two GF(4) Weyl operators commute iff this is 0.  a,b,ap,bp are length-n GF(4) vectors.
    """
    a, b, ap, bp = GF4(a), GF4(b), GF4(ap), GF4(bp)
    s = GF4(0)
    for i in range(len(a)):
        s = s + a[i] * bp[i] + ap[i] * b[i]
    return trace(s)


# --------------------------------------------------------------------------- n-qudit operators
def pauli_tensor(a, b) -> np.ndarray:
    """(X_{a_i} Z_{b_i})_i as a dense 4^n x 4^n matrix."""
    n = len(a)
    op = np.array([[1.0 + 0j]])
    for i in range(n):
        op = np.kron(op, weyl(a[i], b[i]))
    return op


def _self_phase(a, b) -> int:
    """sum_i Tr(a_i b_i): the operator g=(X Z)^{tensor} squares to (-1)^{this} I."""
    a, b = GF4(a), GF4(b)
    return sum(trace(a[i] * b[i]) for i in range(len(a))) % 2


def hermitian_stab(a, b) -> np.ndarray:
    """Hermitian order-2 stabilizer operator for row (a|b).

    g0 = tensor(X_{a_i} Z_{b_i}) satisfies g0^2 = (-1)^{sum Tr(a_i b_i)} I.  Multiply by i^{that}
    to get a Hermitian involution g with g^2 = I and g^dag = g (so (I+g)/2 is a projector).
    """
    g0 = pauli_tensor(a, b)
    s = _self_phase(a, b)
    g = (1j ** s) * g0
    return g


# --------------------------------------------------------------------------- codespace
def codespace_projector(rows, n: int) -> np.ndarray:
    """Projector onto the simultaneous +1 eigenspace of Hermitian stabilizers for symplectic `rows`.

    rows: list of (a,b) with a,b length-n GF(4) vectors, pairwise trace-symplectic-orthogonal.
    Built as the product of commuting projectors (I+g)/2.
    """
    dim = 4 ** n
    Pi = np.eye(dim, dtype=complex)
    for (a, b) in rows:
        g = hermitian_stab(a, b)
        Pi = Pi @ ((np.eye(dim, dtype=complex) + g) / 2.0)
    return Pi


def codespace_basis(rows, n: int):
    """Orthonormal basis (columns) of the codespace from the projector's +1 eigenvectors."""
    Pi = codespace_projector(rows, n)
    H = (Pi + Pi.conj().T) / 2
    w, V = np.linalg.eigh(H)
    cols = V[:, w > 0.5]
    return cols, Pi


def transversal_test(rows, G: np.ndarray, n: int):
    """(leakage, logical, kdim) for transversal G^{tensor n} on the GF(4) code."""
    cols, _ = codespace_basis(rows, n)
    kdim = cols.shape[1]
    if kdim == 0:
        return None, None, 0
    Gn = np.array([[1.0 + 0j]])
    for _ in range(n):
        Gn = np.kron(Gn, G)
    GP = Gn @ cols
    leak = np.linalg.norm(GP - cols @ (cols.conj().T @ GP))
    logical = cols.conj().T @ GP
    return float(leak), logical, kdim


# --------------------------------------------------------------------------- candidate magic gate
def cs_gate() -> np.ndarray:
    """The candidate GF(4) single-qudit level-3 diagonal gate diag(1,1,1,i)."""
    return np.diag([1.0, 1.0, 1.0, 1j]).astype(complex)
