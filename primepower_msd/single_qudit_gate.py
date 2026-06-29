"""Single-qudit diagonal phase-gate builders and a magic-certification routine.

A candidate single-qudit non-Clifford resource is a diagonal gate  diag(zeta_N^{f(x)})  with
zeta_N = exp(2*pi*i/N) and f a low-degree integer polynomial of the canonical lift x in {0,..,d-1}.
``certify_magic`` returns its exact Clifford-hierarchy level plus the substrate anti-collapse
certificate, so the M1 census can decide — operationally, not by assumption — whether a genuine
level-3 single-qudit magic gate exists at d = 2^k.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .clifford_ring import is_clifford, is_level3, level_of
from .weyl import pauli_shift_order


def diag_phase(exponents: list[int], N: int) -> np.ndarray:
    """Diagonal gate diag(exp(2*pi*i * exponents[x] / N)) for x = 0 .. len(exponents)-1."""
    return np.diag([np.exp(2j * np.pi * (e % N) / N) for e in exponents])


def monomial_phase_gate(d: int, degree: int, coeff: int, N: int) -> np.ndarray:
    """Single-qudit gate diag(zeta_N^{coeff * x^degree}), x the canonical lift in {0,..,d-1}.

    Examples (d a power of two):
      * phase_S      ~ monomial_phase_gate(d, 2, 1, 2d)   (Clifford)
      * quadratic T  = monomial_phase_gate(d, 2, 1, 4d)   (the literature-map candidate)
      * cubic T      = monomial_phase_gate(d, 3, 1, 2d)   (non-additive cubic candidate)
    """
    return diag_phase([coeff * (x ** degree) for x in range(d)], N)


@dataclass(frozen=True)
class MagicCertificate:
    d: int
    level: int                 # exact Clifford-hierarchy level (1/2/3/>=4)
    is_clifford: bool
    is_strict_level3: bool     # in C_3 but not in C_2  => single-qudit magic for MSD
    pauli_shift_order: int     # = d; anti-collapse witness (must exceed 2 for k >= 2)
    anticollapse: bool         # cyclic Pauli group not isomorphic to a tensor of qubit groups
    note: str = ""

    @property
    def is_magic(self) -> bool:
        """A genuine single-qudit MSD resource: strictly level 3 and irreducible to qubits."""
        return self.is_strict_level3 and self.anticollapse


def certify_magic(U: np.ndarray, d: int, note: str = "") -> MagicCertificate:
    """Certify a single-qudit gate U on the cyclic ring Z_d: hierarchy level + anti-collapse."""
    lvl = level_of(U, d)
    cliff = is_clifford(U, d)
    strict3 = is_level3(U, d) and not cliff
    order = pauli_shift_order(d)
    return MagicCertificate(
        d=d,
        level=lvl,
        is_clifford=cliff,
        is_strict_level3=strict3,
        pauli_shift_order=order,
        anticollapse=order > 2,
        note=note,
    )
