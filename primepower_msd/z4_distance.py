"""Exact minimum distance of Z_d (d = 2^k) stabilizer codes by syndrome-collision meet-in-the-middle.

Problem.  A logical operator of weight w is a vector v in Z_d^n (CSS Z-type / X-type) or in Z_d^{2n}
(non-CSS, symplectic weight) with  H v == 0 (mod d)  (it commutes with the stabilizers) and v NOT in
the stabilizer module.  d(C) = min such w.  Design copied from qmsd/mindist.py (F_p) with the ring
specifics:

  * syndromes mod d are packed EXACTLY into a uint64 key (bits-per-row = log2 d) when
    rows * log2(d) <= 64; beyond that a polynomial hash is used and every collision is re-verified
    exactly, so there are no false positives either way;
  * weight-a left table (a <= 3) over all coefficient tuples in (Z_d \\ 0)^a (CSS) or
    ((Z_d^2) \\ 0)^a (non-CSS), weight-b right stream (b <= 3) with the leading coefficient normalized
    up to the UNITS of Z_d only (2 is not a unit, so units {1, 3} at d = 4 give 2 classes per
    nonzero coefficient pair; this still enumerates every vector up to a unit scalar, and v is a
    logical iff u v is);
  * a collision certifies  H v == 0  only; v may be a STABILIZER.  Each collision is checked against
    the stabilizer module with a Howell SpanTester and scanning CONTINUES on failure (no early exit
    on the first collision — the lightest element of H^perp is often a stabilizer in CSS codes).

Memory guard: `estimate_left_table` refuses tables beyond `max_left_entries` (default 5e7 ~ 1.5 GB).

API
    css_distance(MX, MZ, n, d=4, d_max=6)        -> (d_min, d_X, d_Z, certified)
    symplectic_distance(rows, n, d=4, d_max=6)   -> (d_min, certified)
    min_weight_logical(H, excl_gens, n, kind, d, d_max) -> (w or None, witness)
"""

from __future__ import annotations

import itertools as it
import math

import numpy as np

from .ringlinalg import howell_form
from .symplectic import SpanTester, symplectic_weight


# --------------------------------------------------------------------------- coefficient alphabets
def _alphabet(kind: str, d: int) -> list[tuple[int, ...]]:
    """Nonzero per-position coefficients: CSS -> (c,), non-CSS -> (a_i, b_i) != (0,0)."""
    if kind == "css":
        return [(c,) for c in range(1, d)]
    return [t for t in it.product(range(d), repeat=2) if t != (0, 0)]


def _unit_reps(kind: str, d: int) -> list[tuple[int, ...]]:
    """One representative per orbit of the alphabet under multiplication by units of Z_d."""
    units = [u for u in range(1, d) if math.gcd(u, d) == 1]
    seen, reps = set(), []
    for t in _alphabet(kind, d):
        if t in seen:
            continue
        reps.append(t)
        for u in units:
            seen.add(tuple((u * x) % d for x in t))
    return reps


def _columns(H: np.ndarray, n: int, kind: str, d: int) -> list[np.ndarray]:
    """Per-position column blocks of H: CSS -> (r,1); non-CSS -> (r,2) = [H[:, i], H[:, n+i]]."""
    if kind == "css":
        return [H[:, i:i + 1] for i in range(n)]
    return [np.stack([H[:, i], H[:, n + i]], axis=1) for i in range(n)]


# --------------------------------------------------------------------------- syndrome keys
class _Keyer:
    def __init__(self, r: int, d: int):
        self.bits = d.bit_length() - 1
        self.r = r
        self.exact = r * self.bits <= 63
        if self.exact:
            self.w = (1 << (self.bits * np.arange(r, dtype=np.uint64))).astype(np.uint64)
        else:
            rng = np.random.default_rng(12345)
            self.w = rng.integers(1, 1 << 62, size=r, dtype=np.uint64) | np.uint64(1)

    def keys(self, S: np.ndarray) -> np.ndarray:
        """S: (N, r) syndromes mod d -> (N,) uint64 keys."""
        S = S.astype(np.uint64)
        if self.exact:
            return (S * self.w[None, :]).sum(axis=1, dtype=np.uint64)
        return (S * self.w[None, :]).sum(axis=1, dtype=np.uint64)   # wraparound hash; verified later


# --------------------------------------------------------------------------- tables
def estimate_left_table(n: int, a: int, kind: str, d: int) -> int:
    return math.comb(n, a) * (len(_alphabet(kind, d)) ** a)


def _vectors_of_weight(n: int, w: int, alphabet, first_reps=None):
    """Yield (positions, coeff tuples) for all weight-w vectors; if first_reps is given the coefficient
    at the FIRST (lowest) position is restricted to those representatives."""
    for pos in it.combinations(range(n), w):
        if first_reps is None:
            for cs in it.product(alphabet, repeat=w):
                yield pos, cs
        else:
            for c0 in first_reps:
                for cs in it.product(alphabet, repeat=w - 1):
                    yield pos, (c0,) + cs


def _syndromes_batch(cols, pos_list, coef_list, r, d):
    """Vectorized syndromes for many (positions, coefficients) of equal weight."""
    N = len(pos_list)
    S = np.zeros((N, r), dtype=np.int64)
    P = np.array(pos_list, dtype=np.int64)          # (N, w)
    Cf = np.array(coef_list, dtype=np.int64)        # (N, w, c)
    for t in range(P.shape[1]):
        blocks = np.stack([cols[p] for p in P[:, t]], axis=0)     # (N, r, c)
        S += np.einsum("nrc,nc->nr", blocks, Cf[:, t, :])
    return S % d


def _build_left(cols, n, a, r, d, kind, keyer, chunk=200000):
    alphabet = _alphabet(kind, d)
    keys, supports = [], []
    buf_pos, buf_cs = [], []

    def flush():
        if not buf_pos:
            return
        S = _syndromes_batch(cols, buf_pos, buf_cs, r, d)
        keys.append(keyer.keys(S))
        supports.append((np.array(buf_pos, dtype=np.int32), np.array(buf_cs, dtype=np.int8)))
        buf_pos.clear(); buf_cs.clear()

    if a == 0:
        keys.append(np.zeros(1, dtype=np.uint64))
        supports.append((np.zeros((1, 0), dtype=np.int32), np.zeros((1, 0, len(alphabet[0])), dtype=np.int8)))
    else:
        for pos, cs in _vectors_of_weight(n, a, alphabet):
            buf_pos.append(pos); buf_cs.append(cs)
            if len(buf_pos) >= chunk:
                flush()
        flush()
    K = np.concatenate(keys)
    order = np.argsort(K, kind="stable")
    P = np.concatenate([s[0] for s in supports]) if supports[0][0].shape[1] else np.zeros((len(K), 0), dtype=np.int32)
    C = np.concatenate([s[1] for s in supports])
    return K[order], P[order], C[order]


def _assemble(n, kind, d, posL, csL, posR, csR):
    v = np.zeros(2 * n if kind != "css" else n, dtype=np.int64)
    for pos, cs in ((posL, csL), (posR, csR)):
        for p, c in zip(pos, cs):
            if kind == "css":
                v[p] = (v[p] + int(c[0])) % d
            else:
                v[p] = (v[p] + int(c[0])) % d
                v[n + p] = (v[n + p] + int(c[1])) % d
    return v


def weight_w_logical(H: np.ndarray, excl_gens, n: int, w: int, kind: str = "css", d: int = 4,
                     max_left_entries: int = 50_000_000, chunk: int = 100000):
    """Find a weight-w vector v with H v == 0 (mod d) not in <excl_gens>; returns (v or None)."""
    a = w // 2
    b = w - a
    if estimate_left_table(n, a, kind, d) > max_left_entries:
        raise MemoryError(f"left table {estimate_left_table(n, a, kind, d)} entries exceeds guard")
    H = np.atleast_2d(np.asarray(H, dtype=np.int64)) % d
    r = H.shape[0]
    cols = _columns(H, n, kind, d)
    keyer = _Keyer(r, d)
    KL, PL, CL = _build_left(cols, n, a, r, d, kind, keyer)
    tester = SpanTester(excl_gens, d) if len(excl_gens) else None
    alphabet = _alphabet(kind, d)
    reps = _unit_reps(kind, d)
    width = len(alphabet[0])

    def check_batch(posR, csR):
        S = _syndromes_batch(cols, posR, csR, r, d)
        KR = keyer.keys((-S) % d)                        # need syndrome(left) == -syndrome(right)
        lo = np.searchsorted(KL, KR, side="left")
        hi = np.searchsorted(KL, KR, side="right")
        for idx in np.nonzero(hi > lo)[0]:
            pr = posR[idx]
            for j in range(lo[idx], hi[idx]):
                pl = PL[j]
                if len(set(pl) & set(pr)):             # positions must be disjoint for weight w
                    continue
                v = _assemble(n, kind, d, pl, CL[j], pr, csR[idx])
                if np.any((H @ v) % d):               # hash collision (non-exact keyer) -> skip
                    continue
                if tester is None or not tester.contains(v):
                    return v
        return None

    buf_pos, buf_cs = [], []
    for pos, cs in _vectors_of_weight(n, b, alphabet, first_reps=reps):
        buf_pos.append(pos); buf_cs.append(cs)
        if len(buf_pos) >= chunk:
            v = check_batch(buf_pos, buf_cs)
            if v is not None:
                return v
            buf_pos, buf_cs = [], []
    if buf_pos:
        return check_batch(buf_pos, buf_cs)
    return None


def min_weight_logical(H, excl_gens, n: int, kind: str = "css", d: int = 4, d_max: int = 6, **kw):
    """(w, witness) for the lightest logical with weight <= d_max, else (None, None)."""
    for w in range(1, d_max + 1):
        v = weight_w_logical(H, excl_gens, n, w, kind, d, **kw)
        if v is not None:
            return w, v
    return None, None


# --------------------------------------------------------------------------- public API
def css_distance(MX, MZ, n: int, d: int = 4, d_max: int = 6, sides=("X", "Z"), **kw):
    """CSS Z_d code with X-stabilizer generators MX and Z-stabilizer generators MZ.
    d_Z = min weight of z in MX^perp \\ <MZ>  (Z-type logicals);  d_X symmetric.
    Returns (d_min, d_X, d_Z, certified) with None entries when > d_max."""
    MX = np.atleast_2d(np.asarray(MX, dtype=np.int64)) % d if len(MX) else np.zeros((0, n), np.int64)
    MZ = np.atleast_2d(np.asarray(MZ, dtype=np.int64)) % d if len(MZ) else np.zeros((0, n), np.int64)
    dZ = dX = None
    if "Z" in sides:
        dZ, _ = min_weight_logical(MX if MX.size else np.zeros((1, n), np.int64), MZ, n, "css", d, d_max, **kw)
    if "X" in sides:
        dX, _ = min_weight_logical(MZ if MZ.size else np.zeros((1, n), np.int64), MX, n, "css", d, d_max, **kw)
    cands = [x for x in (dX, dZ) if x is not None]
    dmin = min(cands) if cands else None
    # an unknown COMPUTED side means > d_max >= dmin; a skipped side leaves dmin uncertified as a min
    certified = dmin is not None and set(sides) == {"X", "Z"}
    return dmin, dX, dZ, certified


def symplectic_distance(rows, n: int, d: int = 4, d_max: int = 6, **kw):
    """Non-CSS: min symplectic weight in C(S) \\ S.  H = rows . J so that H v = 0 <=> v commutes."""
    R = np.atleast_2d(np.asarray(rows, dtype=np.int64)) % d
    H = np.concatenate([R[:, n:], (-R[:, :n]) % d], axis=1) % d      # v.(b | -a) = 0
    w, v = min_weight_logical(H, R, n, "noncss", d, d_max, **kw)
    return w, w is not None


def code_distance(code, d_max: int = 6, sides=("X", "Z"), **kw):
    """Dispatch on a symplectic.Z4Code: (d_min, certified, details).  sides=("Z",) computes only the
    Z-logical distance (useful when M_Z is huge and the X-logicals are known to be heavy)."""
    if code.is_css:
        MX = np.array(code.x_gens, dtype=np.int64).reshape(-1, code.n)
        MZ = np.array(code.z_rows, dtype=np.int64).reshape(-1, code.n)
        dmin, dX, dZ, cert = css_distance(MX, MZ, code.n, code.d, d_max, sides=sides, **kw)
        return dmin, cert, {"d_X": dX, "d_Z": dZ}
    w, cert = symplectic_distance(code.rows, code.n, code.d, d_max, **kw)
    return w, cert, {}
