"""Diagonal gates on m ring qudits Z_d (d = 2^k) as PHASE TABLES, with an exact finite-difference
Clifford-hierarchy oracle, interaction (1/2/3-body) decomposition and a level-3 enumerator.

Convention.  A diagonal gate on m qudits is  D = diag(zeta_M^{P(x)}),  P : Z_d^m -> Z_M,  M = 4d
(d=4: 16th roots, matching the Phase 0-5 carrier q = [0,1,4,13], hierarchy_search's signatures and
ring_transversal.extract_phase).  Embeddings into the table language:
    Pauli   Z^b          ->  P(x) = 4 b x            (mod M)
    Clifford S^a Z^b     ->  P(x) = 2 (a x^2 + 2 b x)  (mod M), evaluated on the lifts x in {0..d-1}
    mod-4 phase i^{F(x)} ->  P(x) = 4 F(x)           (d = 4)
A table is a legitimate gate whatever its values, so the "translation-consistency" problem of
polynomial formulas mod 8/16 on Z_4 inputs disappears: Delta_a P(x) := P(x+a) - P(x) with x+a taken
mod d is always well defined.

Finite-difference hierarchy lemma (d = 2^k, M = 4d; unit vectors e_j):
    P is Pauli-diagonal     <=>  every Delta_{e_j} P is constant;
    P is Clifford-diagonal  <=>  every Delta_{e_i} Delta_{e_j} P is constant;
    P is level-3            <=>  every Delta_{e_i} Delta_{e_j} Delta_{e_k} P is constant.
Proof sketch.  D X_j D^dag = diag(zeta^{Delta_{e_j}P(x - e_j)}) X_j, so D is Clifford iff each
Delta_{e_j}P is a Pauli-diagonal table, and level-3 iff each Delta_{e_j}P is Clifford-diagonal.
A table with all unit differences constant is  P(x) = sum_j c_j x_j + c_0  as a table, and wrapping
x_j -> x_j + d forces d c_j = 0 mod 4d, i.e. c_j in 4Z: that is exactly Z^{c_j/4}.  Conversely
Pauli tables have constant unit differences.  Induct on the level.  General a in Z_d^m need not be
checked because Delta_a is a Z-combination of translates of unit differences and translates of a
constant are constant.  (Cross-validated against the dense oracles in tests/test_phase_tables.py.)
"""

from __future__ import annotations

import itertools as it
from typing import Iterable

import numpy as np

from .ringlinalg import enumerate_module, module_size, right_kernel_howell

Table = np.ndarray


# --------------------------------------------------------------------------- basics
def modulus(d: int) -> int:
    return 4 * d


def zeros(m: int, d: int = 4) -> Table:
    return np.zeros((d,) * m, dtype=np.int64)


def normalize(P: Table, d: int = 4) -> Table:
    """Subtract P(0) (global phase) and reduce mod M."""
    P = np.asarray(P, dtype=np.int64)
    return (P - P[(0,) * P.ndim]) % modulus(d)


def shift(P: Table, a) -> Table:
    """The table x -> P(x + a) (x + a taken mod d coordinatewise)."""
    P = np.asarray(P)
    out = P
    for axis, s in enumerate(a):
        if s % P.shape[axis]:
            out = np.roll(out, -int(s), axis=axis)
    return out


def delta(P: Table, a, d: int = 4) -> Table:
    """Delta_a P (x) = P(x + a) - P(x)  mod M."""
    return (shift(P, a) - np.asarray(P)) % modulus(d)


def iterated_delta(P: Table, alist: Iterable, d: int = 4) -> Table:
    out = np.asarray(P, dtype=np.int64)
    for a in alist:
        out = delta(out, a, d)
    return out


def unit(m: int, j: int, s: int = 1):
    e = [0] * m
    e[j] = s
    return e


# --------------------------------------------------------------------------- hierarchy level
def _all_unit_diffs_constant(P: Table, order: int, d: int) -> bool:
    m = P.ndim
    for axes in it.combinations_with_replacement(range(m), order):
        D = iterated_delta(P, [unit(m, j) for j in axes], d)
        if np.any(D != D.flat[0]):
            return False
    return True


def is_constant(P: Table) -> bool:
    P = np.asarray(P)
    return bool(np.all(P == P.flat[0]))


def is_pauli_table(P: Table, d: int = 4) -> bool:
    return _all_unit_diffs_constant(np.asarray(P, dtype=np.int64), 1, d)


def is_clifford_table(P: Table, d: int = 4) -> bool:
    return _all_unit_diffs_constant(np.asarray(P, dtype=np.int64), 2, d)


def is_level3_table(P: Table, d: int = 4) -> bool:
    return _all_unit_diffs_constant(np.asarray(P, dtype=np.int64), 3, d)


def level(P: Table, d: int = 4) -> int:
    """Clifford-hierarchy level of diag(zeta_M^P): 0 identity (global phase), 1 Pauli, 2 Clifford,
    3 third level, 4 meaning '>= 4'."""
    P = np.asarray(P, dtype=np.int64) % modulus(d)
    if is_constant(P):
        return 0
    for lv in (1, 2, 3):
        if _all_unit_diffs_constant(P, lv, d):
            return lv
    return 4


# --------------------------------------------------------------------------- constructors
def pauli_table(b, d: int = 4) -> Table:
    """Z^{b} on m qudits: P(x) = 4 sum_j b_j x_j."""
    m = len(b)
    grid = np.indices((d,) * m)
    return (4 * sum(int(bj) * grid[j] for j, bj in enumerate(b))) % modulus(d)


def clifford_table(a: int, b: int, d: int = 4) -> Table:
    """Single-qudit S^a Z^b = diag(zeta_{2d}^{a x^2 + 2 b x}) as a table mod 4d."""
    x = np.arange(d)
    return (2 * (a * x * x + 2 * b * x)) % modulus(d)


def from_qudit_tables(tables, d: int = 4) -> Table:
    """Sum_j q_j(x_j) on m qudits (None = identity on that qudit)."""
    m = len(tables)
    out = zeros(m, d)
    for j, q in enumerate(tables):
        if q is None:
            continue
        q = np.asarray(q, dtype=np.int64) % modulus(d)
        shape = [1] * m; shape[j] = d
        out = (out + q.reshape(shape)) % modulus(d)
    return out


def from_poly_mod4(coeffs, mons, m: int) -> Table:
    """diag(i^{F(x)}), F = sum coeffs[t] * prod_{idx in mons[t]} x_idx  (mod 4), as a table mod 16 (d=4)."""
    grid = np.indices((4,) * m)
    F = np.zeros((4,) * m, dtype=np.int64)
    for c, mon in zip(coeffs, mons):
        if int(c) % 4 == 0:
            continue
        term = np.ones((4,) * m, dtype=np.int64)
        for idx in mon:
            term = term * grid[idx]
        F = F + int(c) * term
    return (4 * (F % 4)) % 16


def cz_layer(m: int, pairs, d: int = 4) -> Table:
    """Ring-CS layer: i^{sum c x_i x_j}  (d=4 convention: mod-4 phase -> 4*(.) mod 16).  pairs = [(i,j,c)]."""
    grid = np.indices((d,) * m)
    F = np.zeros((d,) * m, dtype=np.int64)
    for (i, j, c) in pairs:
        F = F + int(c) * grid[i] * grid[j]
    return (4 * (F % d)) % modulus(d)


def ccz_layer(m: int, triples, d: int = 4) -> Table:
    """Ring-CCZ layer: i^{sum c x_i x_j x_k}.  triples = [(i,j,k,c)]."""
    grid = np.indices((d,) * m)
    F = np.zeros((d,) * m, dtype=np.int64)
    for (i, j, k, c) in triples:
        F = F + int(c) * grid[i] * grid[j] * grid[k]
    return (4 * (F % d)) % modulus(d)


def to_unitary(P: Table, d: int = 4) -> np.ndarray:
    """Dense diag(zeta_M^P) with the basis index  idx = sum_j x_j d^{m-1-j}  (big-endian, as the
    dexp_gate / noncss_transversal tensor convention)."""
    P = np.asarray(P, dtype=np.int64) % modulus(d)
    return np.diag(np.exp(2j * np.pi * P.reshape(-1) / modulus(d)))


def from_unitary(U: np.ndarray, m: int, d: int = 4) -> Table:
    """Inverse of to_unitary for a diagonal U whose phases are M-th roots (up to global phase)."""
    M = modulus(d)
    ang = np.angle(np.diag(U)) / (2 * np.pi)
    exps = np.array([int(round(a * M)) % M for a in ang], dtype=np.int64)
    if not np.allclose(np.diag(U), np.exp(2j * np.pi * exps / M), atol=1e-7):
        raise ValueError("U is not diagonal with M-th-root phases")
    return normalize(exps.reshape((d,) * m), d)


# --------------------------------------------------------------------------- dense oracle (validation only)
def _hs_pauli_basis(m: int, d: int):
    from .noncss_transversal import pauli_tensor
    return [pauli_tensor(a, b, d) / np.sqrt(d ** m)
            for a in it.product(range(d), repeat=m) for b in it.product(range(d), repeat=m)]


def level_of_dense(U: np.ndarray, m: int, d: int = 4, tol: float = 1e-8) -> int:
    """Dense multi-qudit hierarchy level via Hilbert-Schmidt Pauli expansion (d^m <= 64).
    U is Pauli iff its expansion has one term; Clifford iff it conjugates the generators X_i, Z_i to
    Paulis; level-3 iff U X_i U^dag X_i^dag and U Z_i U^dag Z_i^dag are Clifford."""
    from .noncss_transversal import pauli_tensor
    basis = _hs_pauli_basis(m, d)
    gens = []
    for i in range(m):
        a = [0] * m; a[i] = 1
        gens.append(pauli_tensor(a, [0] * m, d)); gens.append(pauli_tensor([0] * m, a, d))

    def is_pauli(M):
        c = np.array([np.vdot(Q, M) for Q in basis])
        return int(np.sum(np.abs(c) > tol)) == 1

    def is_cliff(V):
        return all(is_pauli(V @ g @ V.conj().T) for g in gens)

    I = np.eye(U.shape[0])
    if np.allclose(U / U[0, 0], I, atol=tol) and abs(abs(U[0, 0]) - 1) < tol:
        return 0
    if is_pauli(U):
        return 1
    if is_cliff(U):
        return 2
    if all(is_cliff(U @ g @ U.conj().T @ g.conj().T) for g in gens):
        return 3
    return 4


# --------------------------------------------------------------------------- the strict level-3 single-qudit gates
def strict_level3_single(d: int = 4) -> list[tuple[int, ...]]:
    """ALL strict level-3 single-qudit tables with P(0)=0 (d <= 8 exhaustive: M^(d-1) tables)."""
    M = modulus(d)
    out = []
    for rest in it.product(range(M), repeat=d - 1):
        P = np.array((0,) + rest, dtype=np.int64)
        if is_level3_table(P, d) and not is_clifford_table(P, d):
            out.append(tuple(int(v) for v in P))
    return out


def _carriers_d4():
    """The 12 distinct strict level-3 tables produced by hierarchy_search (antidifferences of S^a Z^b),
    recomputed here with the table oracle; [0,1,4,13] is among them."""
    from .hierarchy_search import search
    from .ring_transversal import extract_phase
    out = []
    for a, b, sig, U in search(4):
        phi, N = extract_phase(U, 4)
        assert N in (8, 16)
        tab = tuple((v * (16 // N)) % 16 for v in phi)
        if tab not in out:
            out.append(tab)
    return out


try:  # computed once at import (12 tiny dense checks)
    CARRIERS_D4: list[tuple[int, int, int, int]] = _carriers_d4()
except Exception:  # pragma: no cover - keeps import robust if dense deps move
    CARRIERS_D4 = [(0, 1, 4, 13)]
CARRIER_Q = (0, 1, 4, 13)
assert CARRIER_Q in CARRIERS_D4 or CARRIERS_D4 == [CARRIER_Q]


# --------------------------------------------------------------------------- interaction decomposition
def restrict_to(P: Table, T: Iterable[int]) -> Table:
    """The table x -> P(x_T): coordinates outside T set to 0 (broadcast back to the full grid)."""
    P = np.asarray(P)
    m = P.ndim
    T = set(T)
    idx = tuple(slice(None) if j in T else slice(0, 1) for j in range(m))
    sub = P[idx]
    return np.broadcast_to(sub, P.shape).copy()


def interaction_decomposition(P: Table, d: int = 4) -> dict[frozenset, Table]:
    """Zero-anchored Moebius decomposition  P(x) = sum_{S subset [m]} P_S(x),
    P_S(x) = sum_{T subset S} (-1)^{|S|-|T|} P(x_T);  P_S depends only on x_S and vanishes when any
    coordinate in S is 0.  Returns {frozenset(S): P_S} for all S (including the empty set = P(0))."""
    P = np.asarray(P, dtype=np.int64) % modulus(d)
    m = P.ndim
    M = modulus(d)
    R = {}
    for r in range(m + 1):
        for T in it.combinations(range(m), r):
            R[frozenset(T)] = restrict_to(P, T)
    out = {}
    for S in R:
        acc = np.zeros_like(P)
        for T in R:
            if T <= S:
                sign = -1 if (len(S) - len(T)) % 2 else 1
                acc = acc + sign * R[T]
        out[S] = acc % M
    return out


def one_body(P: Table, d: int = 4) -> Table:
    dec = interaction_decomposition(P, d)
    acc = np.zeros_like(np.asarray(P))
    for S, PS in dec.items():
        if len(S) == 1:
            acc = acc + PS
    return acc % modulus(d)


def mixed_part(P: Table, d: int = 4) -> Table:
    dec = interaction_decomposition(P, d)
    acc = np.zeros_like(np.asarray(P))
    for S, PS in dec.items():
        if len(S) >= 2:
            acc = acc + PS
    return acc % modulus(d)


def is_entangling(P: Table, d: int = 4) -> bool:
    """Not a product of single-qudit diagonals (up to global phase): some |S|>=2 part is nonzero."""
    return bool(np.any(mixed_part(P, d)))


def qudit_marginals(P: Table, d: int = 4) -> list[Table]:
    """The 1-body tables q_j(x_j) = P(x_j e_j) - P(0)."""
    P = np.asarray(P, dtype=np.int64)
    m = P.ndim
    out = []
    for j in range(m):
        idx = [0] * m
        idx[j] = slice(None)
        out.append((P[tuple(idx)] - P[(0,) * m]) % modulus(d))
    return out


def odd_third_difference(P: Table, d: int = 4) -> bool:
    """Some unit third difference is a UNIT multiple of d (= 4 or 12 mod 16 at d=4): the gate has
    genuine order-d (not 2-torsion) cubic content (genuine_v3 criterion, basis independent)."""
    P = np.asarray(P, dtype=np.int64)
    m = P.ndim
    for axes in it.combinations_with_replacement(range(m), 3):
        D = iterated_delta(P, [unit(m, j) for j in axes], d) % modulus(d)
        vals = set(int(v) for v in D.flat)
        for v in vals:
            if v % d == 0 and (v // d) % 2 == 1:
                return True
    return False


def classify(P: Table, d: int = 4) -> dict:
    """Full classification of a diagonal gate table (see module docstring)."""
    from .a2 import escapes_A2, escapes_all_sectors   # local import: a2 depends on this module

    P = normalize(P, d)
    m = P.ndim
    M = modulus(d)
    dec = interaction_decomposition(P, d)
    marg = qudit_marginals(P, d)
    one = [{"qudit": j, "table": [int(v) for v in q], "level": level(q, d),
            "escapes": escapes_all_sectors(q, d)} for j, q in enumerate(marg)]
    grid = np.indices((d,) * m)
    cs, ccz = [], []
    mixed_mod4 = True
    for S, PS in dec.items():
        if len(S) < 2 or not np.any(PS):
            continue
        if np.any(PS % d):
            mixed_mod4 = False
        idx = sorted(S)
        if len(S) == 2:
            i, j = idx
            c = int(iterated_delta(PS, [unit(m, i), unit(m, j)], d)[(0,) * m]) // d % d
            pure = bool(np.all(PS == (d * c * grid[i] * grid[j]) % M))
            cs.append({"qudits": (i, j), "coeff": c, "pure_bilinear": pure})
        elif len(S) == 3:
            i, j, k = idx
            c = int(iterated_delta(PS, [unit(m, i), unit(m, j), unit(m, k)], d)[(0,) * m]) // d % d
            pure = bool(np.all(PS == (d * c * grid[i] * grid[j] * grid[k]) % M))
            ccz.append({"qudits": (i, j, k), "coeff": c, "pure_trilinear": pure})
        else:
            cs.append({"qudits": tuple(idx), "coeff": None, "pure_bilinear": False})
    lv = level(P, d)
    esc = escapes_all_sectors(P, d)
    return {
        "m": m, "d": d, "level": lv,
        "entangling": any(len(S) >= 2 and np.any(PS) for S, PS in dec.items()),
        "escapes_A2": escapes_A2(P, d) if d == 4 else esc,
        "escapes_all_sectors": esc,
        "one_body": one,
        "one_body_escapes": any(o["escapes"] for o in one),
        "mixed_part_escapes": escapes_all_sectors(mixed_part(P, d), d),
        "cs_terms": cs, "ccz_terms": ccz,
        "mixed_parts_mod_d": mixed_mod4,
        "odd_third_difference": odd_third_difference(P, d),
        "genuine_T": lv == 3 and esc,
        "genuine_order_d_cubic": lv == 3 and odd_third_difference(P, d),
    }


# --------------------------------------------------------------------------- enumeration of level-3 tables
def _third_diff_rows(m: int, d: int) -> np.ndarray:
    """Linear functionals (over the d^m table entries) for Delta_i Delta_j Delta_k P(x) - (same at 0),
    all axis triples with repetition and all x != 0."""
    M = modulus(d)
    npts = d ** m
    pts = list(it.product(range(d), repeat=m))
    pos = {p: i for i, p in enumerate(pts)}

    def functional(axes, x):
        row = np.zeros(npts, dtype=np.int64)
        for mask in range(8):
            y = list(x); bits = 0
            for t, ax in enumerate(axes):
                if mask >> t & 1:
                    y[ax] = (y[ax] + 1) % d; bits += 1
            row[pos[tuple(y)]] += (-1) ** (3 - bits)
        return row

    rows = []
    for axes in it.combinations_with_replacement(range(m), 3):
        f0 = functional(axes, (0,) * m)
        for x in pts[1:]:
            rows.append((functional(axes, x) - f0) % M)
    return np.array(rows, dtype=np.int64)


def level3_module(m: int, d: int = 4) -> list[list[int]]:
    """Generators (over Z_M) of the module of level-(<=3) tables on m qudits (as flat vectors of
    length d^m, including the constant direction)."""
    A = _third_diff_rows(m, d)
    return [list(map(int, v)) for v in right_kernel_howell(A, modulus(d))]


def level2_module(m: int, d: int = 4) -> list[list[int]]:
    """Generators of the module of Clifford-diagonal (level <= 2) tables on m qudits."""
    M = modulus(d)
    npts = d ** m
    pts = list(it.product(range(d), repeat=m))
    pos = {p: i for i, p in enumerate(pts)}

    def functional(axes, x):
        row = np.zeros(npts, dtype=np.int64)
        for mask in range(4):
            y = list(x); bits = 0
            for t, ax in enumerate(axes):
                if mask >> t & 1:
                    y[ax] = (y[ax] + 1) % d; bits += 1
            row[pos[tuple(y)]] += (-1) ** (2 - bits)
        return row

    rows = []
    for axes in it.combinations_with_replacement(range(m), 2):
        f0 = functional(axes, (0,) * m)
        for x in pts[1:]:
            rows.append((functional(axes, x) - f0) % M)
    return [list(map(int, v)) for v in right_kernel_howell(np.array(rows, dtype=np.int64), M)]


def rank1_module(m: int, d: int = 4) -> list[list[int]]:
    """Generators of W_m = span{ t(l . x) : l in Z_d^m, t a level<=3 single-qudit table } + Clifford
    tables: the diagonal gates that are products of Clifford-frame-rotated SINGLE-qudit level-3 gates
    and Cliffords (exactly synthesizable from Clifford + single-qudit magic, like qubit CS = T T T^dag)."""
    M = modulus(d)
    pts = np.array(list(it.product(range(d), repeat=m)), dtype=np.int64)
    single = np.array(level3_module(1, d), dtype=np.int64)
    gens = []
    for l in it.product(range(d), repeat=m):
        if not any(l):
            continue
        lx = (pts @ np.array(l)) % d
        for g in single:
            gens.append(list(map(int, g[lx])))
    gens += level2_module(m, d)
    from .ringlinalg import howell_form
    H = howell_form(np.array(gens, dtype=np.int64) % M, M)
    return [list(map(int, r)) for r in H if np.any(np.asarray(r) % M)]


def level3_module_size(m: int, d: int = 4) -> int:
    return module_size(level3_module(m, d), modulus(d))


def howell_enumerate(gens, M: int, max_size: int = 1 << 22) -> np.ndarray:
    """All elements of the Z_M-module <gens> as an (size, ncols) array, via the Howell form: every
    element is uniquely sum_i c_i h_i with 0 <= c_i < M / lead(h_i)  (no 16^r brute product)."""
    from .ringlinalg import howell_form
    H = howell_form(np.array(gens, dtype=np.int64), M)
    rows, ranges = [], []
    for h in H:
        nz = np.nonzero(np.asarray(h) % M)[0]
        if len(nz) == 0:
            continue
        rows.append(np.asarray(h, dtype=np.int64) % M)
        ranges.append(M // int(h[nz[0]]))
    size = int(np.prod(ranges)) if ranges else 1
    if size > max_size:
        raise ValueError(f"module too large to enumerate ({size})")
    out = np.zeros((1, len(gens[0])), dtype=np.int64)
    for h, r in zip(rows, ranges):
        mult = (np.arange(r)[:, None] * h[None, :]) % M            # (r, ncols)
        out = (out[:, None, :] + mult[None, :, :]).reshape(-1, out.shape[1]) % M
    return out


def enumerate_level3(m: int, d: int = 4, max_size: int = 1 << 22) -> np.ndarray:
    """All level-(<=3) tables on m qudits with P(0)=0, as a (count, d^m) flat array (exact)."""
    gens = level3_module(m, d)
    E = howell_enumerate(gens, modulus(d), max_size=max_size)
    E = (E - E[:, :1]) % modulus(d)                      # kill the global-phase direction
    E = np.unique(E, axis=0)
    return E


def _shift_index(m: int, d: int, a) -> np.ndarray:
    """Permutation idx' such that (flat table)[idx'] = P(x + a) under big-endian flattening."""
    pts = np.array(list(it.product(range(d), repeat=m)), dtype=np.int64)
    shifted = (pts + np.asarray(a)) % d
    return (shifted * (d ** np.arange(m - 1, -1, -1))).sum(axis=1)


def batch_flags(E: np.ndarray, m: int, d: int = 4) -> dict[str, np.ndarray]:
    """Vectorized classification of many flat tables (rows of E, P(0)=0 not required).
    Returns boolean arrays: level3 (<=3), clifford (<=2), entangling, escapes (all sectors),
    one_body_escapes, mixed_escapes, mixed_mod_d."""
    M = modulus(d)
    E = np.asarray(E, dtype=np.int64) % M
    N = E.shape[0]

    def Dl(T, a):
        return (T[:, _shift_index(m, d, a)] - T) % M

    def all_const(T):
        return np.all(T == T[:, :1], axis=1)

    lvl3 = np.ones(N, bool); cliff = np.ones(N, bool)
    for axes in it.combinations_with_replacement(range(m), 2):
        T = E
        for j in axes:
            T = Dl(T, unit(m, j))
        cliff &= all_const(T)
    for axes in it.combinations_with_replacement(range(m), 3):
        T = E
        for j in axes:
            T = Dl(T, unit(m, j))
        lvl3 &= all_const(T)

    def escapes(T):
        k = d.bit_length() - 1
        esc_all = np.ones(N, bool)
        for j in range(1, k):
            s1, s2 = 2 ** j, d // 2 ** j
            pres = np.ones(N, bool)
            for i in range(m):
                for l in range(m):
                    D = Dl(Dl(T, unit(m, i, s1)), unit(m, l, s2))
                    pres &= np.all(D == 0, axis=1)
            esc_all &= ~pres
        return esc_all

    # one-body part: sum_j (P(x_j e_j) - P(0)); mixed = P - one_body - P(0)
    pts = np.array(list(it.product(range(d), repeat=m)), dtype=np.int64)
    one = np.zeros_like(E)
    for j in range(m):
        proj = np.zeros_like(pts); proj[:, j] = pts[:, j]
        idx = (proj * (d ** np.arange(m - 1, -1, -1))).sum(axis=1)
        one = (one + E[:, idx] - E[:, :1]) % M
    mixed = (E - E[:, :1] - one) % M
    return {"level3": lvl3, "clifford": cliff, "entangling": np.any(mixed != 0, axis=1),
            "escapes": escapes(E), "one_body_escapes": escapes(one), "mixed_escapes": escapes(mixed),
            "mixed_mod_d": np.all(mixed % d == 0, axis=1)}


def sample_level3(m: int, n_samples: int, rng: np.random.Generator, d: int = 4) -> list[Table]:
    gens = np.array(level3_module(m, d), dtype=np.int64)
    M = modulus(d)
    out = []
    for _ in range(n_samples):
        c = rng.integers(0, M, size=len(gens))
        P = (c @ gens) % M
        out.append(normalize(P.reshape((d,) * m), d))
    return out


def catalog(m: int, d: int = 4, n_samples: int = 0, seed: int = 0) -> dict:
    """E1: classify every (or a sample of) level-3 table on m qudits; check the 1-body lemma
    'escapes => the 1-body part escapes' and 'mixed parts are multiples of d'."""
    rng = np.random.default_rng(seed)
    size = level3_module_size(m, d)
    if n_samples == 0 and size <= (1 << 22):
        E = enumerate_level3(m, d)
        mode = "exhaustive"
    else:
        E = np.array([P.reshape(-1) for P in sample_level3(m, n_samples or 100000, rng, d)])
        mode = "sampled"
    counts = {"total": int(E.shape[0])}
    violations = []
    for start in range(0, E.shape[0], 65536):
        F = batch_flags(E[start:start + 65536], m, d)
        strict = F["level3"] & ~F["clifford"]
        counts["level3_strict"] = counts.get("level3_strict", 0) + int(strict.sum())
        counts["entangling"] = counts.get("entangling", 0) + int((strict & F["entangling"]).sum())
        esc = strict & F["escapes"]
        counts["escaping"] = counts.get("escaping", 0) + int(esc.sum())
        counts["escaping_entangling"] = counts.get("escaping_entangling", 0) + int((esc & F["entangling"]).sum())
        viol = esc & (~F["one_body_escapes"] | F["mixed_escapes"])
        counts["lemma_violations"] = counts.get("lemma_violations", 0) + int(viol.sum())
        counts["mixed_not_mod_d"] = counts.get("mixed_not_mod_d", 0) + int((strict & ~F["mixed_mod_d"]).sum())
        for row in E[start:start + 65536][viol][:10 - len(violations)]:
            violations.append([int(v) for v in row])
    return {"m": m, "d": d, "mode": mode, "module_size": int(size), "counts": counts,
            "lemma_violation_examples": violations}
