"""Existence search: are there genuine single-qudit level-3 (C_3 \\ C_2) gates on the cyclic
ring Z_{2^k}?  (Kill-criterion K1, settled constructively rather than by a naive monomial scan.)

The M1 census shows the *naive* monomial phase gates diag(zeta_N^{x^j}) are NOT level 3 for
d = 2^k (k >= 2) — the cyclic wraparound defect pushes them to level >= 4.  This is a sharper
(corrected) statement than the literature map.  But it does not settle existence.

Constructive handle.  For a diagonal U with cyclic phase antidifference  U X U^dag = X . diag(v),
one has U in C_3  iff  diag(v) in C_2.  So: take ANY diagonal Clifford V = S^a Z^b, rescale its
global phase so the phases close around the cycle, integrate (antidifference) to U, and test
whether U is strictly level 3.  The diagonal Clifford group is finite ({S^a Z^b}), so this is an
exact, exhaustive existence test that scales to d = 32.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from primepower_msd.weyl import clock_Z, clock_X, Zpow, omega
from primepower_msd.clifford_ring import phase_S, is_clifford, is_level3, level_of
from primepower_msd.single_qudit_gate import certify_magic

DIMS = [2, 4, 8, 16, 32]


def diagonal_clifford(d: int, a: int, b: int) -> np.ndarray:
    """Diagonal Clifford S^a Z^b = diag(zeta_{2d}^{a x^2 + 2 b x})."""
    z2d = np.exp(2j * np.pi / (2 * d))
    return np.diag([z2d ** ((a * x * x + 2 * b * x) % (2 * d)) for x in range(d)])


def antidifference(V: np.ndarray, d: int) -> np.ndarray:
    """Integrate a diagonal V around the cycle to U with  U X U^dag = X . V'  (V' = phase-rescaled V).

    Rescale the global phase of V so the product of its diagonal is 1 (necessary for the cyclic
    phases to close), then set p_0 = 1, p_{x+1} = p_x * v_x.  Returns U = diag(p).
    """
    v = np.diag(V).astype(complex)
    prod = np.prod(v)
    lam = prod ** (-1.0 / d)          # principal d-th root; makes prod(lam*v) = 1
    v = lam * v
    p = np.ones(d, dtype=complex)
    for x in range(d - 1):
        p[x + 1] = p[x] * v[x]
    return np.diag(p)


def search(d: int):
    """Return list of strict level-3 single-qudit diagonal gates found at dimension d (as phase tuples)."""
    found = []
    seen = set()
    for a in range(2 * d):
        for b in range(d):
            V = diagonal_clifford(d, a, b)
            U = antidifference(V, d)
            if is_level3(U, d) and not is_clifford(U, d):
                # canonical signature: phase exponents relative to p_0, rounded
                phases = np.angle(np.diag(U)) / (2 * np.pi)
                sig = tuple(int(round((ph - phases[0]) * (4 * d))) % (4 * d) for ph in phases)
                if sig not in seen:
                    seen.add(sig)
                    found.append((a, b, sig, U))
    return found


def main() -> int:
    print("=" * 78)
    print("K1 EXISTENCE SEARCH — single-qudit C_3 \\ C_2 gates on cyclic Z_{2^k}")
    print("  (antidifference of diagonal Cliffords; exact + exhaustive over S^a Z^b)")
    print("=" * 78)
    verdict = {}
    for d in DIMS:
        hits = search(d)
        verdict[d] = len(hits)
        k = d.bit_length() - 1
        print(f"\n  d = {d:<3} (Z_2^{k}): {len(hits)} distinct strict-level-3 single-qudit gate(s) found")
        if hits:
            a, b, sig, U = hits[0]
            cert = certify_magic(U, d, note=f"antidiff(S^{a} Z^{b})")
            diagvals = np.diag(U)
            shown = ", ".join(f"{v.real:+.3f}{v.imag:+.3f}i" for v in diagvals[: min(d, 6)])
            print(f"        example: antidiff(S^{a} Z^{b})  ->  level {cert.level}, "
                  f"magic={cert.is_magic}, anticollapse={cert.anticollapse}")
            print(f"        phase exps (units 2pi/{4*d}): {sig}")
            print(f"        diag[:6]: {shown}{' ...' if d > 6 else ''}")

    print("\n" + "=" * 78)
    print("VERDICT (K1)")
    print("=" * 78)
    for d in DIMS:
        if d == 2:
            continue
        if verdict[d] > 0:
            print(f"  d = {d:<3}:  single-qudit level-3 magic gate EXISTS ({verdict[d]} found) "
                  f"-> K1 CLEARED")
        else:
            print(f"  d = {d:<3}:  NO single-qudit C_3 gate -> K1 FIRES (pivot / no-go)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
