"""primepower_msd — single-qudit magic state distillation over the cyclic ring Z_{2^k}.

Phase 0 substrate (parametric in d = 2^k, k >= 1): the cyclic clock Heisenberg-Weyl
group, its single-qudit Clifford group, an exact Clifford-hierarchy level oracle, and a
single-qudit magic-gate census. Everything is built from exact d x d complex matrices with
root-of-unity phases, so even-d phase-convention subtleties are handled by direct computation
rather than by hand.

Target program: see RESEARCH_PLAN.md / LITERATURE_MAP.md in this folder.
"""

from .weyl import omega, clock_X, clock_Z, pauli, pauli_shift_order
from .clifford_ring import (
    phase_S,
    fourier_H,
    multiplier_M,
    is_pauli,
    is_clifford,
    is_level3,
    level_of,
)
from .single_qudit_gate import diag_phase, monomial_phase_gate, certify_magic

__all__ = [
    "omega",
    "clock_X",
    "clock_Z",
    "pauli",
    "pauli_shift_order",
    "phase_S",
    "fourier_H",
    "multiplier_M",
    "is_pauli",
    "is_clifford",
    "is_level3",
    "level_of",
    "diag_phase",
    "monomial_phase_gate",
    "certify_magic",
]
