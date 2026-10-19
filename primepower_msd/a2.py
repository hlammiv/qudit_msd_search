"""Universality criterion for ring qudits Z_d, d = 2^k: sector (A_j) escape.

Background (Borda-Rincon-Galindo, arXiv:2512.20787).  For d = p^m the Clifford group's adjoint action
on sl(d^n) is REDUCIBLE: the Pauli operators X^u Z^v split into invariant sectors W_j indexed by
gcd(u, v, p^m) = p^j.  Clifford + G is universal only if G couples the sectors (orbit mixing).
For d = 2^k, j = 1..k-1, the sector algebras

    A_j := span{ X^a Z^b : all a_i, b_i in 2^j Z_d }          (A_1 = "A2" at d = 4)

are the obstructions: A_1 at d = 4 is a maximal ABELIAN subalgebra (X^2, Z^2 commute), i.e. one fixed
product basis {(|0>+-|2>), (|1>+-|3>)}, and its normalizer N(A_1) = {monomial unitaries in that basis}
is a closed proper subgroup of U(4^n).  Hence a gate that preserves every A_j under conjugation can
never make Clifford + G dense.

Diagonal criterion (exact).  For D = diag(zeta_M^P), M = 4d, and a shift by u = 2^j e_i,
    D X^u D^dag = diag(zeta_M^{Delta_u P(x - u)}) X^u.
This lies in A_j iff the diagonal factor is a combination of Z^{2^j v}, i.e. a function of
x mod d/2^j, i.e. invariant under shifts by (d/2^j) e_l for every l:
    D normalizes A_j  <=>  Delta_{2^j e_i} Delta_{(d/2^j) e_l} P == 0 (mod M)  for all i, l.
(Z^{2^j v}-type diagonals commute with every Pauli power in A_j, so the Z-side is automatic.)
At d = 4 this is  q(x + 2) - q(x) in {0, 8} for all x  (A2 preserved), the criterion stated in
UNIVERSALITY_CHECK.md.  Dense cross-checks are in tests/test_a2.py.

Known values (d = 4): Cliffords, the ring-CCZ i^{xyz}, i^{x^3} and every diag(i^F) with F a
polynomial mod 4 PRESERVE A2 (non-universalizing); the carrier zeta16^[0,1,4,13], zeta8^{x^3},
zeta16^{x} ESCAPE.
"""

from __future__ import annotations

import itertools as it

import numpy as np

from .phase_tables import iterated_delta, modulus, normalize, to_unitary, unit


# --------------------------------------------------------------------------- symbolic (diagonal tables)
def sector_defect(P, j: int, d: int = 4) -> np.ndarray:
    """Max over (i,l) of |Delta_{2^j e_i} Delta_{(d/2^j) e_l} P| as a table of residues (0 = preserved)."""
    P = np.asarray(P, dtype=np.int64)
    m = P.ndim
    M = modulus(d)
    s1, s2 = 2 ** j, d // (2 ** j)
    worst = np.zeros_like(P)
    for i in range(m):
        for l in range(m):
            D = iterated_delta(P, [unit(m, i, s1), unit(m, l, s2)], d) % M
            worst = np.maximum(worst, D)
    return worst


def preserves_sector(P, j: int = 1, d: int = 4) -> bool:
    """True iff diag(zeta_M^P) normalizes A_j (1 <= j <= k-1)."""
    k = d.bit_length() - 1
    assert 1 <= j <= k - 1, "sectors exist only for 1 <= j <= k-1"
    return not np.any(sector_defect(P, j, d))


def escapes_A2(P, d: int = 4) -> bool:
    """d = 4 shorthand: does the gate leave A2 = A_1?"""
    return not preserves_sector(P, 1, d)


def escapes_all_sectors(P, d: int = 4) -> bool:
    """Universality-enabling (diagonal) iff every sector A_j, j = 1..k-1, is escaped."""
    k = d.bit_length() - 1
    if k <= 1:
        return True            # qubits: no proper sector; any non-Clifford suffices
    return all(not preserves_sector(P, j, d) for j in range(1, k))


def escape_witness(P, d: int = 4):
    """(x, j, i, l, residue) for the first violated sector condition, or None if all preserved."""
    P = np.asarray(P, dtype=np.int64)
    m = P.ndim
    M = modulus(d)
    k = d.bit_length() - 1
    for j in range(1, k):
        s1, s2 = 2 ** j, d // (2 ** j)
        for i in range(m):
            for l in range(m):
                D = iterated_delta(P, [unit(m, i, s1), unit(m, l, s2)], d) % M
                nz = np.argwhere(D != 0)
                if len(nz):
                    x = tuple(int(t) for t in nz[0])
                    return {"x": x, "j": j, "i": i, "l": l, "residue": int(D[x])}
    return None


# --------------------------------------------------------------------------- dense (any small unitary)
def _pauli_basis(m: int, d: int, step: int = 1):
    """Dense X^a Z^b for all a,b in (step*Z_d)^m, normalized to unit Hilbert-Schmidt norm."""
    from .noncss_transversal import pauli_tensor
    rng = range(0, d, step)
    ops = []
    for a in it.product(rng, repeat=m):
        for b in it.product(rng, repeat=m):
            ops.append(pauli_tensor(a, b, d) / np.sqrt(d ** m))
    return ops


def preserves_sector_dense(U: np.ndarray, m: int, d: int = 4, j: int = 1, tol: float = 1e-9) -> bool:
    """Dense check that U A_j U^dag subset A_j (generators X_i^{2^j}, Z_i^{2^j} suffice since
    conjugation is an algebra automorphism).  Feasible for d^m <= 64."""
    from .noncss_transversal import pauli_tensor
    k = d.bit_length() - 1
    assert 1 <= j <= k - 1
    full = _pauli_basis(m, d, 1)
    sector_keys = set()
    rng = range(0, d, 2 ** j)
    for a in it.product(rng, repeat=m):
        for b in it.product(rng, repeat=m):
            sector_keys.add((a, b))
    outside = []
    idx = 0
    for a in it.product(range(d), repeat=m):
        for b in it.product(range(d), repeat=m):
            if (a, b) not in sector_keys:
                outside.append(full[idx])
            idx += 1
    Ud = U.conj().T
    gens = []
    for i in range(m):
        a = [0] * m; a[i] = 2 ** j
        gens.append(pauli_tensor(a, [0] * m, d))
        gens.append(pauli_tensor([0] * m, a, d))
    for g in gens:
        C = U @ g @ Ud
        for Q in outside:
            if abs(np.vdot(Q, C)) > tol:
                return False
    return True


def escapes_A2_dense(U: np.ndarray, m: int, d: int = 4, tol: float = 1e-9) -> bool:
    return not preserves_sector_dense(U, m, d, 1, tol)


def escapes_all_sectors_dense(U: np.ndarray, m: int, d: int = 4, tol: float = 1e-9) -> bool:
    k = d.bit_length() - 1
    return all(not preserves_sector_dense(U, m, d, j, tol) for j in range(1, k))


def table_escapes_dense(P, d: int = 4) -> bool:
    P = normalize(P, d)
    return escapes_all_sectors_dense(to_unitary(P, d), P.ndim, d)
