"""Distance certificate for ring-punctured Reed–Muller codes (Z-side).

Setting: P ⊂ F_2^m, code = lift of RM(r,m) punctured at P (k = |P| when phi_r(P) is independent).
A Z-logical of weight w  <=>  S = A ∪ R in RM(r,m)^perp with A ⊆ P nonempty, R ∩ P = ∅, |R| = w
(Z_4 words reduce to binary: a residue in the binary stabilizer has weight >= 2^{r+1}).

Lemma (half-space reduction).  If A ⊄ aff(R), pick an affine form l with l|_R = 0 and l = 1 somewhere on A.
For every polynomial g of degree <= r-1, sum_S l*g = sum_{A ∩ {l=1}} g = 0, so A ∩ {l=1} is a nonempty
RM(r-1,m)^perp word inside P ∩ {l=1}.  Hence if
    (H)  phi_{r-1}(P ∩ H) is linearly independent for every affine half-space H = {l = 1},
every logical has A ⊆ aff(R) =: F, dim F <= |R| - 1, and S is a word of RM(r, dim F)^perp on F.
So d_Z >= D follows from (H) plus a local statement on flats of dimension <= D-2 (case 2), which for
Sidon-type P is settled by Sidon bounds s(k) and the minimum-weight / Kasami–Tokura structure of
RM(dim F - r - 1, dim F).

This module checks (H) exactly, verifies Sidon/half-Sidon properties, ranks, and builds explicit
flat witnesses for the matching upper bound.
"""

from __future__ import annotations

import itertools as it

import numpy as np


def monomials(m: int, deg: int):
    return [S for k in range(deg + 1) for S in it.combinations(range(m), k)]


def phi_matrix(pts: np.ndarray, deg: int) -> np.ndarray:
    """rows = points, cols = monomials of degree <= deg (uint8 0/1)."""
    pts = np.asarray(pts, dtype=np.uint8)
    cols = []
    for S in monomials(pts.shape[1], deg):
        c = np.ones(len(pts), dtype=np.uint8)
        for i in S:
            c &= pts[:, i]
        cols.append(c)
    return np.stack(cols, axis=1) if cols else np.zeros((len(pts), 0), np.uint8)


def gf2_rank_rows(M: np.ndarray) -> int:
    """GF(2) rank, rows packed into Python ints."""
    rows = [int("".join(map(str, r)), 2) for r in np.asarray(M, dtype=np.uint8)]
    basis = {}
    rank = 0
    for v in rows:
        while v:
            h = v.bit_length() - 1
            if h in basis:
                v ^= basis[h]
            else:
                basis[h] = v
                rank += 1
                break
    return rank


def to_ints(pts) -> list[int]:
    return [int(sum(int(b) << i for i, b in enumerate(p))) for p in pts]


def is_sidon(ints: list[int]) -> bool:
    seen = set()
    for a, b in it.combinations(ints, 2):
        s = a ^ b
        if s in seen:
            return False
        seen.add(s)
    return True


def halfspace_check(pts: np.ndarray, deg: int) -> dict:
    """(H): rank phi_deg(P ∩ {a.x + b = 1}) == |P ∩ ...| for all a != 0, b in {0,1}."""
    pts = np.asarray(pts, dtype=np.uint8)
    m = pts.shape[1]
    Phi = phi_matrix(pts, deg)
    ints = np.array(to_ints(pts), dtype=np.int64)
    worst = None
    max_size = 0
    fails = 0
    for a in range(1, 2 ** m):
        dot = np.array([bin(int(x) & a).count("1") & 1 for x in ints])
        for b in (0, 1):
            sel = (dot ^ b) == 1
            sz = int(sel.sum())
            max_size = max(max_size, sz)
            rk = gf2_rank_rows(Phi[sel]) if sz else 0
            if rk < sz:
                fails += 1
                if worst is None:
                    worst = {"a": a, "b": b, "size": sz, "rank": rk}
    return {"deg": deg, "halfspaces": 2 * (2 ** m - 1), "max_|P∩H|": max_size, "fails": fails,
            "first_fail": worst, "ok": fails == 0}


def affine_span(points_int: list[int]) -> list[int]:
    base = points_int[0]
    dirs = []
    piv = {}
    for p in points_int[1:]:
        v = p ^ base
        while v:
            h = v.bit_length() - 1
            if h in piv:
                v ^= piv[h]
            else:
                piv[h] = v
                dirs.append(v)
                break
    span = [base]
    for d in dirs:
        span += [x ^ d for x in span]
    return span


def flat_witness(P_int: list[int], subset: list[int]) -> dict:
    """The affine span F of the given P-points; returns |F|, |F ∩ P| and R = F \\ P."""
    F = affine_span(subset)
    Ps = set(P_int)
    R = [x for x in F if x not in Ps]
    return {"dim": int(np.log2(len(F))), "size": len(F), "occupancy": len(F) - len(R), "R": R}


# --------------------------------------------------------------------------- refined (flat-wise) certificate
def left_kernel(M: np.ndarray) -> np.ndarray:
    """Basis (rows) of {x : x M = 0} over GF(2)."""
    M = np.array(M, dtype=np.uint8) % 2
    n, c = M.shape
    A = np.concatenate([M, np.eye(n, dtype=np.uint8)], axis=1)
    r = 0
    for col in range(c):
        piv = np.nonzero(A[r:, col])[0]
        if len(piv) == 0:
            continue
        p = r + piv[0]
        if p != r:
            A[[r, p]] = A[[p, r]]
        rows = np.nonzero(A[:, col])[0]
        rows = rows[rows != r]
        A[rows] ^= A[r]
        r += 1
        if r == n:
            break
    return A[r:, c:]


def form_values(ints: np.ndarray, form: int, m: int) -> np.ndarray:
    """affine form encoded as a | (c << m): value a.x + c on each point."""
    a, c = form & ((1 << m) - 1), form >> m
    v = np.zeros(len(ints), dtype=np.uint8)
    x = ints & a
    while np.any(x):
        v ^= (x & 1).astype(np.uint8)
        x >>= 1
    return v ^ c


def dependent_points_by_halfspace(P: np.ndarray, deg: int = 2) -> dict[int, np.ndarray]:
    """For every nonconstant affine form l: boolean mask of the points of P ∩ {l = 1} that lie in the
    support of some phi_deg-dependency inside P ∩ {l = 1} (only forms with a nonempty mask are kept)."""
    P = np.asarray(P, dtype=np.uint8)
    m = P.shape[1]
    ints = np.array(to_ints(P), dtype=np.int64)
    Phi = phi_matrix(P, deg)
    out = {}
    for a in range(1, 2 ** m):
        for c in (0, 1):
            form = a | (c << m)
            sel = np.nonzero(form_values(ints, form, m))[0]
            if len(sel) == 0:
                continue
            K = left_kernel(Phi[sel])
            if len(K) == 0:
                continue
            mask = np.zeros(len(P), dtype=bool)
            mask[sel[np.any(K, axis=0)]] = True
            out[form] = mask
    return out


def affine_cosets_in(Fset, k: int, nbits: int, limit: int = 10 ** 6):
    """Distinct cosets l0 + W (W linear of dim k) contained in Fset (a set of nbits-bit ints).
    Bitmask DFS: cand[x] = 1 iff x + span ⊆ Fset; adding w updates cand <- cand & cand[x ^ w].
    Yields sorted tuples of the 2^k coset elements."""
    N = 1 << nbits
    T = np.zeros(N, dtype=bool)
    T[list(Fset)] = True
    idx = np.arange(N)
    seen = set()

    def rec(l0, span, cand, last):
        if len(span) == 2 ** k:
            key = tuple(sorted(l0 ^ s for s in span))
            if key not in seen:
                seen.add(key)
                yield key
            return
        for x in np.nonzero(cand)[0]:
            w = int(x) ^ l0
            if w <= last or w in span:
                continue
            c2 = cand & cand[idx ^ w]
            if not c2[l0]:
                continue
            yield from rec(l0, span | {s ^ w for s in span}, c2, w)
            if len(seen) >= limit:
                return

    for l0 in sorted(Fset):
        yield from rec(l0, {0}, T.copy(), 0)
        if len(seen) >= limit:
            return


def flat_from_forms(forms: list[int], m: int) -> list[int] | None:
    """Common zero set of affine forms (as point ints); None if empty."""
    pts = np.arange(2 ** m, dtype=np.int64)
    ok = np.ones(2 ** m, dtype=bool)
    for f in forms:
        ok &= form_values(pts, f, m) == 0
    z = np.nonzero(ok)[0]
    return [int(x) for x in z] if len(z) else None


def exact_flat_check(P: np.ndarray, forms: list[int], r: int = 3) -> dict:
    """G = {all forms = 0}.  Returns the dimension of X_G = {B ⊆ P \\ G : B ⊥ (cubics vanishing on G)}
    (0  <=>  every cubic dependency S with S \\ P ⊆ G has S ⊆ G) and a basis of X_G."""
    P = np.asarray(P, dtype=np.uint8)
    m = P.shape[1]
    ints = np.array(to_ints(P), dtype=np.int64)
    vals = np.stack([form_values(ints, f, m) for f in forms], axis=1)
    off = np.nonzero(vals.any(axis=1))[0]
    Phi = phi_matrix(P[off], r - 1)
    cols = [Phi * vals[off, i:i + 1] for i in range(len(forms))]
    M = np.concatenate(cols, axis=1)
    K = left_kernel(M)
    basis = []
    for row in K:
        v = np.zeros(len(P), dtype=np.uint8)
        v[off] = row
        basis.append(v)
    return {"off_G": int(len(off)), "dim_X": int(len(K)), "basis": basis}


# --------------------------------------------------------------------------- top-down failing-flat BFS
class FlatBFS:
    """Exact enumeration of the flats G (given by a space W_G of affine forms) for which
    X_G = {B ⊆ P \\ G : B ⊥ I_G} is nonzero, I_G = cubics vanishing on G = sum_{l in W_G} l * RM(r-1).

    X_G = 0  <=>  rank{ val_l & mono : l in basis(W_G), mono of degree <= r-1 } == |P \\ G|
    (rows are bitmasks over P).  If G (codim c) fails then some codim-(c-1) flat containing it fails
    (B' ⊄ G'' for some of the 2^c - 1 flats G'' ⊇ G of codim c-1), so the failing flats are found
    top-down: children of failing flats only."""

    def __init__(self, P: np.ndarray, r: int = 3):
        self.P = np.asarray(P, dtype=np.uint8)
        self.n, self.m = self.P.shape
        ints = np.array(to_ints(self.P), dtype=np.int64)
        self.mono = [int("".join(map(str, col[::-1])), 2) for col in phi_matrix(self.P, r - 1).T]
        self.val = {}
        for f in range(1, 2 ** (self.m + 1)):
            if f == (1 << self.m):
                continue
            v = form_values(ints, f, self.m)
            self.val[f] = int("".join(map(str, v[::-1])), 2)

    @staticmethod
    def _insert(basis: dict, v: int) -> bool:
        while v:
            h = v.bit_length() - 1
            b = basis.get(h)
            if b is None:
                basis[h] = v
                return True
            v ^= b
        return False

    def rows_for(self, f: int):
        vf = self.val[f]
        return [vf & mu for mu in self.mono]

    def check(self, forms: list[int], parent_basis: dict | None = None):
        """Return (fails, basis, outside) for W = span(forms); reuses a parent's row basis for forms[:-1]."""
        basis = dict(parent_basis) if parent_basis is not None else {}
        todo = forms[-1:] if parent_basis is not None else forms
        for f in todo:
            for row in self.rows_for(f):
                self._insert(basis, row)
        outside = 0
        for f in forms:
            outside |= self.val[f]
        return len(basis) < bin(outside).count("1"), basis, outside


def canon_forms(forms: list[int]) -> tuple:
    """Reduced echelon basis of span(forms) (affine forms as ints) — canonical key."""
    basis = {}
    for v in forms:
        while v:
            h = v.bit_length() - 1
            if h in basis:
                v ^= basis[h]
            else:
                basis[h] = v
                break
    hs = sorted(basis, reverse=True)
    for h in hs:
        for h2 in hs:
            if h2 != h and (basis[h2] >> h) & 1:
                basis[h2] ^= basis[h]
    return tuple(sorted(basis.values()))
