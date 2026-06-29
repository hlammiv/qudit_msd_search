"""Phase 0 driver: M0 substrate calibration + M1 single-qudit magic-gate census.

Run:  python -m primepower_msd.census       (from the prime_msd repo root)
  or:  python primepower_msd/census.py

M0 validates the Clifford-hierarchy oracle against the textbook qubit (d=2) answers and the
three verified algebraic facts.  M1 then scans, for d = 2, 4, 8, 16, 32, a family of candidate
single-qudit diagonal phase gates and reports the EXACT Clifford-hierarchy level of each, plus the
anti-collapse certificate — answering operationally whether a genuine level-3 single-qudit magic
gate exists at each power-of-two dimension (kill-criterion K1).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if __package__ in (None, ""):  # allow `python primepower_msd/census.py`
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from primepower_msd.clifford_ring import phase_S, fourier_H, is_clifford, level_of
from primepower_msd.ring import squaring_is_additive, cubing_is_additive
from primepower_msd.single_qudit_gate import monomial_phase_gate, diag_phase, certify_magic

DIMS = [2, 4, 8, 16, 32]
_LVL = {1: "Pauli", 2: "Clifford", 3: "LEVEL-3", 4: ">=4/none"}


def gauss_sum(d: int) -> complex:
    """Quadratic Gauss sum G(d) = sum_x exp(2*pi*i x^2 / d)."""
    return sum(np.exp(2j * np.pi * (x * x) / d) for x in range(d))


def m0_calibration() -> bool:
    """Validate the oracle and the verified algebraic facts. Returns True iff all pass."""
    print("=" * 78)
    print("M0 — SUBSTRATE CALIBRATION")
    print("=" * 78)
    ok = True

    # (1) Qubit ground truth: S Clifford, T level-3, H Clifford.
    s2 = phase_S(2)                                   # diag(1, i)
    t2 = monomial_phase_gate(2, 2, 1, 4 * 2)          # diag(1, e^{i pi/4}) = qubit T
    h2 = fourier_H(2)
    checks = [
        ("d=2  S = diag(1,i)        is Clifford (lvl 2)", level_of(s2, 2) == 2),
        ("d=2  T = diag(1,e^{ipi/4}) is LEVEL-3 (lvl 3)", level_of(t2, 2) == 3),
        ("d=2  H (Hadamard)          is Clifford (lvl 2)", level_of(h2, 2) == 2),
    ]

    # (2) diag(zeta_{2d}^{x^2}) is the Clifford S for every d (the corrected precision fact).
    for d in DIMS:
        checks.append((f"d={d:<2} diag(zeta_{2*d}^(x^2)) is Clifford (= S, NOT magic)",
                       is_clifford(phase_S(d), d)))

    # (3) Non-additivity of squaring/cubing over Z_{2^k}, k>=2 (the anti-collapse lever).
    for d in DIMS:
        if d == 2:
            continue
        checks.append((f"d={d:<2} x->x^2 NON-additive over Z_{d}", not squaring_is_additive(d)))
        checks.append((f"d={d:<2} x->x^3 NON-additive over Z_{d}", not cubing_is_additive(d)))

    # (4) Quadratic Gauss sum G(d) = (1+i) sqrt(d) for d = 0 mod 4.
    for d in DIMS:
        if d % 4 == 0:
            G = gauss_sum(d)
            checks.append((f"d={d:<2} Gauss sum sum_x e^(2pi i x^2/{d}) = (1+i)sqrt({d})",
                           np.isclose(G, (1 + 1j) * np.sqrt(d), atol=1e-9)))

    for label, passed in checks:
        print(f"  [{'PASS' if passed else 'FAIL'}]  {label}")
        ok = ok and passed
    print(f"\n  M0 {'PASSED — oracle trustworthy.' if ok else 'FAILED — STOP, fix the oracle.'}\n")
    return ok


def candidate_gates(d: int):
    """The candidate single-qudit diagonal gates scanned at dimension d.

    Yields (name, U, expectation_note).
    """
    yield ("S = quad @ 2d   diag(z_{2d}^{x^2})", phase_S(d), "calibration: Clifford")
    yield ("quad @ 4d       diag(z_{4d}^{x^2})", monomial_phase_gate(d, 2, 1, 4 * d), "lit-map candidate")
    yield ("quad @ 8d       diag(z_{8d}^{x^2})", monomial_phase_gate(d, 2, 1, 8 * d), "higher precision")
    yield ("cubic @ 2d      diag(z_{2d}^{x^3})", monomial_phase_gate(d, 3, 1, 2 * d), "non-additive cubic")
    yield ("cubic @ 4d      diag(z_{4d}^{x^3})", monomial_phase_gate(d, 3, 1, 4 * d), "non-additive cubic")
    yield ("linear @ 2d     diag(z_{2d}^{x})  ", monomial_phase_gate(d, 1, 1, 2 * d), "Z^{1/2}")
    # BRG-style linear universality gate T_s = diag(e^{2pi i x / s}); non-Clifford iff s does not
    # divide K_d = 2d. Take s = 4d (4d never divides 2d for d>=1).
    yield ("BRG T_s (s=4d)  diag(e^{2pi i x/4d})", diag_phase(list(range(d)), 4 * d), "universality target")


def m1_census() -> dict[int, bool]:
    """Census the Clifford-hierarchy level of candidate single-qudit gates. Returns {d: has_magic}."""
    print("=" * 78)
    print("M1 — SINGLE-QUDIT MAGIC-GATE CENSUS  (level / clifford / strict-L3 / magic)")
    print("=" * 78)
    found = {}
    for d in DIMS:
        k = d.bit_length() - 1
        print(f"\n  d = {d}  (Z_2^{k}),  Pauli shift order = {d}  "
              f"[anti-collapse {'OK' if d > 2 else '(qubit baseline)'}]")
        print(f"    {'candidate':<34} {'level':<10} {'magic?':<7} note")
        print(f"    {'-'*34} {'-'*10} {'-'*7} {'-'*22}")
        any_magic = False
        for name, U, note in candidate_gates(d):
            cert = certify_magic(U, d)
            any_magic = any_magic or cert.is_magic
            print(f"    {name:<34} {_LVL[cert.level]:<10} "
                  f"{('YES' if cert.is_magic else '-'):<7} {note}")
        found[d] = any_magic
    return found


def main() -> int:
    np.set_printoptions(precision=4, suppress=True)
    ok = m0_calibration()
    if not ok:
        print("ABORT: M0 calibration failed; downstream results are not trustworthy.")
        return 1
    found = m1_census()

    print("\n" + "=" * 78)
    print("VERDICT  (K1: does a genuine single-qudit level-3 magic gate exist?)")
    print("=" * 78)
    for d in DIMS:
        if d == 2:
            continue
        verdict = "EXISTS — proceed" if found[d] else "NONE FOUND — K1 risk; widen family or pivot"
        print(f"  d = {d:<3}:  {verdict}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
