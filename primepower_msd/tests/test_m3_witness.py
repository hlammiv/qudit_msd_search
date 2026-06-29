"""Regression test for THE result: a verified single-ququart (d=4) magic-state-distillation code.

A [[9,1,>=2]]_4 ring CSS code over the cyclic clock Z_4 admits a transversal genuinely-quadratic
level-3 gate (antidiff(S^3 Z^1)) inducing the strict level-3 logical gate diag(1,-i,1,i).  Verified
algebraically AND by dense codespace simulation in C^(4^9).

Run:  python -m pytest primepower_msd/tests/test_m3_witness.py -q
"""

import numpy as np

from primepower_msd.ringlinalg import enumerate_module, right_kernel_howell, in_span, module_size
from primepower_msd.clifford_ring import is_level3, is_clifford, level_of
from primepower_msd.ring_transversal import extract_phase
from primepower_msd.hierarchy_search import search as fam_search

D = 4
N_QUDITS = 9
A = np.array([[1, 3, 3, 0, 1, 1, 1, 1, 1]])          # X-stabilizer generator (cyclic Z_4)
G = np.array([1, 0, 0, 0, 1, 1, 2, 0, 3])            # logical-X coset generator
GATE = (3, 1)                                         # transversal physical gate = antidiff(S^3 Z^1)


def _setup():
    fam = fam_search(D)
    gi = next(i for i, (a, b, s, U) in enumerate(fam) if (a, b) == GATE)
    U = fam[gi][3]
    phi, Nmod = extract_phase(U, D)
    phi = list(phi)
    mx = enumerate_module(A, D)
    Lgen = np.vstack([A, G])
    L = enumerate_module(Lgen, D)
    mz = right_kernel_howell(Lgen, D)
    cosets = [tuple((t * G) % D) for t in range(D)]
    return U, phi, Nmod, mx, L, mz, cosets


def _Phi(state, phi, Nmod):
    return sum(phi[int(x)] for x in state) % Nmod


def test_physical_gate_is_genuinely_quadratic():
    """The transversal gate must be one of the HARD (trivial-translation-stabilizer) level-3 gates."""
    _, phi, Nmod, *_ = _setup()
    stab1 = [s for s in range(D)
             if all((phi[(x + s) % D] - phi[x]) % Nmod == (phi[s % D] - phi[0]) % Nmod for x in range(D))]
    assert stab1 == [0]            # genuinely-quadratic, NOT qubit-like


def test_valid_css_k1_cyclic():
    _, phi, Nmod, mx, L, mz, cosets = _setup()
    assert all((A[0] @ np.array(b)) % D == 0 for b in mz)              # X/Z commutation
    assert len(mx) * module_size(np.array(mz), D) == D ** (N_QUDITS - 1)  # k=1
    assert next(t for t in range(1, D + 1) if tuple((t * G) % D) in set(mx)) == D  # cyclic Z_4 logical


def test_distance_at_least_2():
    _, phi, Nmod, mx, L, mz, cosets = _setup()
    mx_set = set(mx)
    assert not any(sum(1 for x in v if x) == 1 for v in L if v not in mx_set)   # no weight-1 X-logical
    z_wt1 = [(i, w) for i in range(N_QUDITS) for w in range(1, D)
             if all((w * int(a[i])) % D == 0 for a in mx)
             and not in_span([0] * i + [w] + [0] * (N_QUDITS - 1 - i), mz, D)]
    assert not z_wt1                                                            # no weight-1 Z-logical


def test_transversal_and_induces_strict_level3():
    U, phi, Nmod, mx, L, mz, cosets = _setup()
    # transversality: total phase constant on each M_X-coset (codewords are |M_X|-sparse, so this is exact)
    assert all(_Phi(tuple((np.array(r) + np.array(a)) % D), phi, Nmod) == _Phi(r, phi, Nmod)
               for r in cosets for a in mx)
    lam0 = _Phi(cosets[0], phi, Nmod)
    lam = [(_Phi(cosets[l], phi, Nmod) - lam0) % Nmod for l in range(D)]
    UL = np.diag([np.exp(2j * np.pi * l / Nmod) for l in lam])
    assert level_of(UL, D) == 3 and not is_clifford(UL, D)         # strict level-3 induced gate


def test_minimal_8_1_2_code():
    """The MINIMAL single-ququart MSD code [[8,1,2]]_4 (the [[9,1,2]] with its free |0> ancilla removed):
    transversal antidiff(S^3 Z^1) induces diag(1,-i,1,i); gamma = log(8)/log(2) = 3.0."""
    import primepower_msd.ringlinalg as RL
    import primepower_msd.ring_distance as RD
    n = 8
    A8 = np.array([[1, 3, 3, 1, 1, 1, 1, 1]])
    g8 = np.array([1, 0, 0, 1, 1, 2, 0, 3])
    fam = fam_search(D)
    U = fam[next(i for i, (a, b, s, U) in enumerate(fam) if (a, b) == GATE)][3]
    phi, Nmod = extract_phase(U, D)
    phi = list(phi)
    mx = RL.enumerate_module(A8, D)
    mx_set = set(mx)
    cosets = [tuple((t * g8) % D) for t in range(D)]
    assert all(_Phi(a, phi, Nmod) == 0 for a in mx)                                   # self
    assert all(_Phi(tuple((np.array(r) + np.array(a)) % D), phi, Nmod) == _Phi(r, phi, Nmod)
               for r in cosets for a in mx)                                           # transversal
    Lg = np.vstack([A8, g8])
    mz = RL.right_kernel_howell(Lg, D)
    dist = RD.min_logical_weight(RL.howell_form(A8, D), mz, n, D, cap=6)
    assert dist == 2
    assert abs(RD.code_gamma(n, 1, dist) - 3.0) < 1e-9
    # dense codespace sim in C^(4^8)
    size = D ** n
    idx = np.arange(size)
    Dvec = np.diag(U).astype(complex)
    Uvec = np.ones(size, dtype=complex)
    for i in range(n):
        Uvec *= Dvec[(idx // (D ** i)) % D]
    enc = lambda v: sum(int(v[i]) * (D ** i) for i in range(n))
    P = np.zeros((D, size), dtype=complex)
    for li, r in enumerate(cosets):
        for a in mx:
            P[li, enc([(r[j] + a[j]) % D for j in range(n)])] += 1.0
        P[li] /= np.linalg.norm(P[li])
    UP = Uvec[None, :] * P
    Lmat = P.conj() @ UP.T
    assert np.linalg.norm(UP.T - P.T @ Lmat) < 1e-9                                   # transversal (leakage 0)
    assert is_level3(Lmat, D) and not is_clifford(Lmat, D)                            # logical strict level-3


def test_dense_codespace_simulation():
    """Gold standard: build the codespace in C^(4^9) and confirm D^{⊗9} preserves it (leakage 0)."""
    U, phi, Nmod, mx, L, mz, cosets = _setup()
    size = D ** N_QUDITS
    idx = np.arange(size)
    Dvec = np.diag(U).astype(complex)
    Uvec = np.ones(size, dtype=complex)
    for i in range(N_QUDITS):
        Uvec *= Dvec[(idx // (D ** i)) % D]
    enc = lambda v: sum(int(v[i]) * (D ** i) for i in range(N_QUDITS))
    P = np.zeros((D, size), dtype=complex)
    for li, r in enumerate(cosets):
        for a in mx:
            P[li, enc([(r[j] + a[j]) % D for j in range(N_QUDITS)])] += 1.0
        P[li] /= np.linalg.norm(P[li])
    UP = Uvec[None, :] * P
    Lmat = P.conj() @ UP.T
    leak = np.linalg.norm(UP.T - P.T @ Lmat)
    assert leak < 1e-9                                              # transversal: no codespace leakage
    assert np.linalg.norm(Lmat - np.diag(np.diag(Lmat))) < 1e-9     # pure diagonal logical gate
    assert is_level3(Lmat, D) and not is_clifford(Lmat, D)          # simulated logical gate is strict level-3
