"""Transversality of a diagonal single-qudit gate on a ring CSS code, and the induced logical gate.

Key reduction (valid because the magic gate is DIAGONAL).  With X-basis logical states
|psi_c> ∝ sum_{a in M_X} |c + a>, the transversal gate U = D^{⊗n}, D = diag(zeta_N^{phi(x)}), gives

    U |psi_c> = sum_{a in M_X} zeta_N^{Phi(c+a)} |c+a>,   Phi(y) = sum_i phi(y_i mod d).

U is diagonal, so it never mixes M_X-cosets; hence U preserves the codespace **iff** the total phase
Phi is constant (mod N) on every M_X-coset of every valid logical rep.  Then the induced logical gate
is the diagonal diag(zeta_N^{Phi(r_ell)}) over the d logical values.  This makes the M2 existence
question a finite phase-arithmetic check — no d^n state vectors required.
"""

from __future__ import annotations

import numpy as np

from .clifford_ring import is_clifford, is_level3, level_of
from .ring_css import CSSCode


def extract_phase(U: np.ndarray, d: int) -> tuple[list[int], int]:
    """Recover (phi, N) of a diagonal gate U = diag(exp(2 pi i phi(x)/N)), phi[0] normalized to 0."""
    ang = np.angle(np.diag(U)) / (2 * np.pi)
    for N in (2 * d, 4 * d, 8 * d, 16 * d):
        exps = [int(round(a * N)) % N for a in ang]
        if np.allclose(np.diag(U), [np.exp(2j * np.pi * (e / N)) for e in exps], atol=1e-7):
            e0 = exps[0]
            return [(e - e0) % N for e in exps], N
    raise ValueError("phase precision exceeds 16d; widen the search")


def stab1(phi: list[int], N: int, d: int) -> dict[int, int]:
    """Single-coordinate translation stabilizer of phi: {s : phi(x+s)-phi(x) is constant in x} -> that constant.

    A subgroup of Z_d (with a homomorphism delta:s->const).  The X-stabilizer translation stabilizer
    Stab(phi) = {a in Z_d^n : each a_i in stab1, sum_i delta(a_i) == 0 mod N} is the set of X-stabilizers
    on which the diagonal gate diag(zeta_N^phi) acts trivially -> globally guarantees transversality.
    """
    out = {}
    for s in range(d):
        deltas = [(phi[(x + s) % d] - phi[x]) % N for x in range(d)]
        if all(v == deltas[0] for v in deltas):
            out[s] = deltas[0]
    return out


def stab_elements(phi: list[int], N: int, d: int, n: int, max_size: int = 1 << 16):
    """Enumerate Stab(phi) subset Z_d^n (X-stabilizer vectors on which the gate acts trivially)."""
    from itertools import product as _product
    s1 = stab1(phi, N, d)
    keys = list(s1.keys())
    if len(keys) ** n > max_size:
        return None
    elems = []
    for a in _product(keys, repeat=n):
        if sum(s1[v] for v in a) % N == 0:
            elems.append(a)
    return elems


def _Phi(state, phi: list[int], N: int, weights=None) -> int:
    if weights is None:
        return sum(phi[int(x)] for x in state) % N
    return sum(w * phi[int(x)] for w, x in zip(weights, state)) % N


def induced_logical(code: CSSCode, phi: list[int], N: int, weights=None) -> np.ndarray | None:
    """Return the induced logical diagonal gate if the (addressable) transversal gate is logical, else None.

    The physical gate is ⊗_i D^{w_i} (uniform if weights is None); requires a k=1 cyclic code and checks
    the coset-constancy (transversality) condition exactly.  Vectorized over the M_X module (the inner
    loop) — Phi(r + a) for all stabilizers a at once — which is the search hot path.
    """
    if code.k != 1 or not code.cyclic:
        return None
    d = code.d
    phi_arr = np.asarray(phi, dtype=np.int64)
    nq = code.n
    w = np.ones(nq, dtype=np.int64) if weights is None else np.asarray(weights, dtype=np.int64)
    mx = np.asarray(code.mx, dtype=np.int64)            # (|M_X|, n)
    lam = []
    for r in code.cosets:
        states = (np.asarray(r, dtype=np.int64)[None, :] + mx) % d      # (|M_X|, n)
        phis = (phi_arr[states] * w[None, :]).sum(axis=1) % N           # (|M_X|,)
        if not np.all(phis == phis[0]):
            return None                                # phase not constant on coset -> not transversal
        lam.append(int(phis[0]))
    lam0 = lam[0]
    return np.diag([np.exp(2j * np.pi * ((l - lam0) % N) / N) for l in lam])


def certify_distillation_code(code: CSSCode, U: np.ndarray, weights=None):
    """Full M2 certificate for (code, physical gate U, optional addressable weights). Verdict dict."""
    d = code.d
    phi, N = extract_phase(U, d)
    UL = induced_logical(code, phi, N, weights)
    if UL is None:
        return {"transversal": False}
    lvl = level_of(UL, d)
    cliff = is_clifford(UL, d)
    strict3 = is_level3(UL, d) and not cliff
    return {
        "transversal": True,
        "logical_level": lvl,
        "logical_is_clifford": cliff,
        "logical_strict_level3": strict3,
        "code_anticollapse": code.cyclic,     # logical qudit is a genuine cyclic Z_d
        "is_distillation_code": bool(strict3 and code.cyclic),
        "logical_phases": list(np.round(np.angle(np.diag(UL)) / (2 * np.pi) * (4 * d)).astype(int) % (4 * d)),
    }
