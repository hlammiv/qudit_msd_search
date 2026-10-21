"""Measurement-based (HHPW-style) distillation over Z_d: inner codes on which the T-check is a transversal Clifford.

Haah–Hastings–Poulin–Wecker (arXiv:1703.07847) distil by *measuring* H = T X T^dag on k encoded noisy magic
states with an inner [[n,k,D]] code (check = T^dag . controlled-X^u . T on the data, 2n noisy T per check) and
arranging the checks with a sensitive outer code.  Lemma 4.4 / Eq. (4.4):  with outer row weight k and
s = (D-1)/2,  noisy T per output  =  D + (D-1)(n/k - 1)  for D-th order suppression, so
        gamma_HHPW = log(D + (D-1)(n/k - 1)) / log D   ->  1  as  n/k -> 1.

Z_d version (cultivation.py): T X T^dag = X S_{2d}^{a} (Pauli) with a odd <=> T escapes, and
T^{-1} = N T N^dag (Clifford) with N: x -> -x (q(-x) + q(x) is Clifford for the escaping carriers), so T and
T^{-1} cost one magic state each.  The transversal check  prod_i T^{tau_i} X^{u_i} T^{-tau_i}  is a logical
Clifford iff the 0/1 X-stabiliser rows S with signs tau in {+-1}^n are *signed weakly self-dual*:
    Q(s) = sum_{i in s} tau_i = 0 (mod 2d),   B(s,t) = sum_{i in s cap t} tau_i = 0 (mod d)      (all s,t in S),
with Z-stabilisers  S diag(tau).  Exact integer balance (= 0 over Z) makes the code work for EVERY d at once.
The logical X module S^perp_tau / S carries the form G = L diag(tau) L^T; a *normal* basis (G = diag(odd)) makes
the check act as prod_j X_j S_j^{odd}: one escaping (T-type) check per logical, i.e. every output is a T state.
Binary reduction: normal <=> some logical has odd weight (HHPW Thm 3.4); then Gram–Schmidt with odd pivots
over Z_d.

Family (exact balance): S = shortened-punctured RM(r,m) (monomials of degree 1..r on F_2^m \\ {0}),
tau(x) = (-1)^{x_1 + ... + x_{2r+1}}  (m >= 2r+1): every product of two generators has degree <= 2r, so its
signed sum over F_2^m vanishes, and the punctured point carries no generator.
"""

from __future__ import annotations

import itertools as it
import math

import numpy as np

from .ringlinalg import right_kernel_howell


def outer_degree(D: int) -> int:
    """checks per input s: (D-1)/2 for D >= 5 (HHPW Lemma 4.15 needs s >= 2, so D = 3 also uses s = 2)."""
    return max(2, (D - 1) // 2)


def t_per_output(n: int, k: int, D: int) -> float:
    """1 input + 2 n T per check, s checks of weight k per input: 1 + 2 s n / k  (= D + (D-1)(n/k-1), D >= 5)."""
    return 1 + 2 * outer_degree(D) * n / k


def gamma_hhpw(n: int, k: int, D: int) -> float | None:
    if k <= 0 or D <= 1:
        return None
    return math.log(t_per_output(n, k, D)) / math.log(D)


# --------------------------------------------------------------------------- families
def signed_rm(r: int, m: int):
    """(S, tau): S rows = monomials x^A, 1 <= |A| <= r, evaluated on F_2^m \\ {0}; tau = (-1)^{x_1+..+x_{2r+1}}."""
    if m < 2 * r + 1:
        raise ValueError("need m >= 2r+1 for exact balance")
    pts = np.array([[(v >> j) & 1 for j in range(m)] for v in range(1, 2 ** m)], dtype=np.int64)
    rows = []
    for deg in range(1, r + 1):
        for A in it.combinations(range(m), deg):
            rows.append(np.prod(pts[:, list(A)], axis=1))
    tau = 1 - 2 * (pts[:, :2 * r + 1].sum(1) % 2)
    return np.array(rows, dtype=np.int64), tau


def balance_report(S, tau, d: int | None = None) -> dict:
    """Max |Q| and |B| over generators (exact integers); with d also the mod-(2d)/(d) verdict."""
    S = np.asarray(S, dtype=np.int64)
    B = (S * tau) @ S.T
    Q = np.diag(B)
    rep = {"max_abs_Q": int(np.abs(Q).max()), "max_abs_B": int(np.abs(B - np.diag(Q)).max()) if len(S) > 1 else 0,
           "exact": bool(not np.any(B))}
    if d is not None:
        rep["ok_mod_d"] = bool(not np.any(Q % (2 * d)) and not np.any(B % d))
    return rep


# --------------------------------------------------------------------------- logicals
def _rank2(A):
    A = np.array(A, dtype=np.int64) % 2
    r = 0
    for c in range(A.shape[1]):
        piv = np.nonzero(A[r:, c])[0]
        if len(piv) == 0:
            continue
        p = r + piv[0]
        A[[r, p]] = A[[p, r]]
        for i in range(A.shape[0]):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        r += 1
        if r == A.shape[0]:
            break
    return r


def logical_x_basis(S, tau, d: int):
    """Z_d vectors L (k x n) in S^perp_tau = {x : S diag(tau) x = 0} completing S (independent mod 2)."""
    S = np.asarray(S, dtype=np.int64) % d
    K = np.array(right_kernel_howell((S * tau) % d, d), dtype=np.int64).reshape(-1, S.shape[1])
    basis2 = [r for r in S % 2]
    rank = _rank2(basis2)
    L = []
    for v in K:
        if not np.any(v % 2):
            continue
        if _rank2(basis2 + [v % 2]) > rank:
            basis2.append(v % 2); rank += 1; L.append(v % d)
    return np.array(L, dtype=np.int64).reshape(-1, S.shape[1])


def _gram(L, tau, d):
    G = (L * tau) @ L.T
    return G % d, np.diag(G) % (2 * d)


def normal_basis(L, tau, d: int):
    """Basis change of the logical X reps to G = diag(odd) mod d (and odd Q mod 2d); None if hyperbolic."""
    L = [np.asarray(v, dtype=np.int64) for v in L]
    ip = lambda a, b: int(((a * tau) @ b))
    done, rest = [], list(L)
    while rest:
        j = next((i for i, v in enumerate(rest) if ip(v, v) % 2), None)
        if j is None:
            if not done:
                return None
            # rest alternating mod 2: I_1 + lambda -> I_3 with e (odd), a, b (a.b odd)
            e = done.pop()
            a = rest.pop(0)
            jb = next((i for i, v in enumerate(rest) if ip(a, v) % 2), None)
            if jb is None:
                return None
            b = rest.pop(jb)
            rest = [e + a, e + b, e + a + b] + rest
            continue
        v = rest.pop(j)
        qv = ip(v, v)
        inv = pow(qv % d, -1, d) if d > 1 else 0
        rest = [(w - (ip(w, v) * inv % d) * v) % d for w in rest]
        done.append(v % d)
    N = np.array(done, dtype=np.int64)
    G, Q = _gram(N, tau, d)
    off = G - np.diag(np.diag(G))
    if np.any(off % d) or not np.all(Q % 2):
        # Gram–Schmidt pass did not converge mod d (should not happen with odd pivots); report failure
        return None
    return N


def analyze(S, tau, d: int = 4, D_max: int = 6, distance: bool = True, label: str = "") -> dict:
    S = np.asarray(S, dtype=np.int64)
    n = S.shape[1]
    r = _rank2(S)
    k = n - 2 * r
    bal = balance_report(S, tau, d)
    out = {"label": label, "n": n, "k": k, "d": d, "balance": bal}
    if not bal["ok_mod_d"]:
        out["status"] = "not signed self-orthogonal"; return out
    L = logical_x_basis(S, tau, d)
    out["k_found"] = int(len(L))
    N = normal_basis(L, tau, d) if len(L) else None
    out["normal"] = N is not None
    if N is not None:
        _, Q = _gram(N, tau, d)
        out["logical_S_exponents"] = sorted(set(int(x) for x in Q % (2 * d)))     # prod S^{tau} acts as S_L^{Q_j}
    if distance:
        from .z4_distance import css_distance
        dmin, dX, dZ, cert = css_distance(S % d, (S * tau) % d, n, d, D_max, sides=("Z",))
        out["d_Z"] = dZ if dZ is not None else f">{D_max}"
        D = dZ if dZ is not None else D_max + 1
        out["gamma_hhpw_bound"] = gamma_hhpw(n, k, D) if k > 0 else None
        out["T_per_output"] = t_per_output(n, k, D) if k > 0 else None
    out["status"] = "ok"
    return out
