"""Cultivation substrate for d = 2^k ring qudits: the T-check is a transversal Clifford on signed 6.6.6
color codes.

Facts (verified in tests/test_cultivation.py):
  * For every strict level-3 diagonal T on Z_d (d = 4, 8, 16), T X T^dag = X diag(zeta_M^g) with g a Clifford
    table c + 2(a x^2 + 2 b x); T escapes all Pauli sectors (universality-enabling)  <=>  a is odd.
    So the cultivation check operator is X S^{odd} (the analogue of H_XY ~ S X for qubits).
  * The signed Z_d triangular color code (X rows 0/1 plaquettes, Z rows plaquettes times the honeycomb
    bipartition tau = +-1) is a CSS code [[ (3D^2+1)/4, 1, D ]]_d; X^{(x)n} is logical X, and the transversal
    product  prod_i T^{tau_i} X T^{-tau_i}  acts as T_L X_L T_L^dag with the SAME carrier T_L (dense check
    on n = 7, d = 4; exact symbolic hosting for D <= 9, d = 4, 8, 16).  Hence |T_L> = T_L|+_L> is a +1
    eigenstate of a transversal physical check: the qubit cultivation check carries over to every 2^k
    (the argument is algebraic: balanced plaquettes and +- overlaps cancel, logical exponent = sum tau = +-1).
"""

from __future__ import annotations

import itertools as it
from math import comb

import numpy as np

from . import phase_tables as PT
from . import rgt
from .a2 import escapes_all_sectors
from .phase_tables import modulus
from .symplectic import Z4Code


def color_code_plaquettes(distance: int):
    """6.6.6 triangular color code from stim's generated circuit: (H, coords) with H the plaquette matrix."""
    import stim
    c = stim.Circuit.generated("color_code:memory_xyz", distance=distance, rounds=2).flattened()
    anc, data, sup = set(), [], {}
    for ins in c:
        t = [x.value for x in ins.targets_copy()]
        if ins.name == "MR":
            anc |= set(t)
        elif ins.name == "MY":
            data = t
        elif ins.name == "CX":
            for a, b in zip(t[::2], t[1::2]):
                sup.setdefault(b, set()).add(a)
    idx = {q: i for i, q in enumerate(data)}
    H = np.zeros((len(anc), len(data)), dtype=np.int64)
    for r, a in enumerate(sorted(anc)):
        for q in sup[a]:
            H[r, idx[q]] = 1
    coords = c.get_final_qubit_coordinates()
    return H, np.array([coords[q] for q in data], dtype=float)


def honeycomb_signs(H: np.ndarray, coords: np.ndarray) -> np.ndarray:
    """Bipartition tau in {+1,-1}^n of the honeycomb (nearest neighbours sharing a plaquette)."""
    n = H.shape[1]
    adj = [[j for j in range(n) if j != i and np.hypot(*(coords[i] - coords[j])) < 1.2 and np.any(H[:, i] & H[:, j])]
           for i in range(n)]
    tau = np.zeros(n, dtype=np.int64)
    tau[0] = 1
    stack = [0]
    while stack:
        i = stack.pop()
        for j in adj[i]:
            if tau[j] == 0:
                tau[j] = -tau[i]
                stack.append(j)
            assert tau[j] == -tau[i], "lattice not bipartite"
    assert np.all(tau) and not np.any(H @ tau), "plaquettes not balanced"
    return tau


def signed_color_code(distance: int, d: int = 4) -> tuple[Z4Code, np.ndarray]:
    H, xy = color_code_plaquettes(distance)
    tau = honeycomb_signs(H, xy)
    MX, MZ = H.copy(), (H * tau) % d
    assert not np.any((MX @ MZ.T) % d)
    return Z4Code.from_css(MX, MZ, H.shape[1], d=d, label=f"CC{distance}_Z{d}"), tau


def escaping_carrier(d: int) -> np.ndarray:
    """Smallest escaping strict level-3 table q = A y + B C(y,2) + C C(y,3) (mod 4d), d-periodic."""
    M = modulus(d)
    for C, B, A in it.product(range(M), repeat=3):
        q = np.array([(A * y + B * comb(y, 2) + C * comb(y, 3)) % M for y in range(2 * d)])
        if np.array_equal(q[:d], q[d:]) and PT.level(q[:d], d) == 3 and escapes_all_sectors(q[:d], d):
            return q[:d]
    raise ValueError(d)


def check_table(q, d: int) -> np.ndarray:
    """g with T X T^dag = X diag(zeta^g) up to the index shift: g(x) = q(x+1) - q(x)."""
    M = modulus(d)
    return np.array([(q[(x + 1) % d] - q[x]) % M for x in range(d)])


def s_exponent(F, d: int) -> int:
    """Quadratic exponent a of a Clifford table F = c + 2(a x^2 + 2 b x) (mod 4d)."""
    M = modulus(d)
    F = np.asarray(F).reshape(-1) % M
    return int(((F[2 % d] - 2 * F[1] + F[0]) % M) // 4)


def check_layer(code: Z4Code, tau, q) -> rgt.TransversalLayer:
    """Diagonal part of prod_i T^{tau_i} X T^{-tau_i}: per-qudit tables tau_i * g."""
    M = modulus(code.d)
    g = check_table(q, code.d)
    return rgt.TransversalLayer.per_qudit([(int(t) * g) % M for t in tau], code.d, "T-check")


def verify_check(distance: int, d: int = 4) -> dict:
    """Exact (symbolic, CSS support mode) verification of the transversal T-check on the color code."""
    code, tau = signed_color_code(distance, d)
    q = escaping_carrier(d)
    layer = check_layer(code, tau, q)
    hosted = bool(rgt.hosts(code, layer, "support").hosted)
    F = rgt.induced_logical(code, layer)
    ones = np.ones(code.n, dtype=np.int64)
    zr = np.array(code.z_rows, dtype=np.int64).reshape(-1, code.n)
    return {"distance": distance, "d": d, "n": code.n, "k": code.k, "hosted": hosted,
            "logical_g": np.asarray(F).reshape(-1).tolist(), "logical_s_exponent": s_exponent(F, d) % (2 * d),
            "physical_g": check_table(q, d).tolist(), "x_all_ones_logical": not np.any((zr @ ones) % d)}
