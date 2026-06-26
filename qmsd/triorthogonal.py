"""Triorthogonal codes and the CSS distillation-code builder (NOTES sec 3, sec 10 item 7).

A linear space C subset F_p^n is *triorthogonal* (NOTES Def 1) iff, for the chosen
basis rows b_a, BOTH the quadratic (self-orthogonality) and cubic conditions hold,
including repeated indices:

    (pair)   sum_i (b_a * b_b)_i       = 0  (mod p)   for all a <= b
    (triple) sum_i (b_a * b_b * b_c)_i = 0  (mod p)   for all a <= b <= c

The CSS magic-state code is built by puncturing a triorthogonal generator G at a set
of k columns (NOTES eq 1, sec 3): the punctured generator G' (= Gp) and the shortened
generator G0 satisfy the nesting G0 subset G' subset G0^perp.  The code is

    CSS(G0 -> X,  G'^perp -> Z)   ==  [[n0 - k, k, d]]_p

with X-stabilizers G0 and Z-stabilizers the dual of the punctured code (NOTES sec 3,
sec 10 item 7).  The [[n0-k, k, d]] label is only valid when the surviving (punctured-
column) submatrix of G is full rank -- always checked.
"""
from __future__ import annotations

import numpy as np

from .field import GFp


def is_triorthogonal(basis, p) -> bool:
    """True iff the F_p span of ``basis`` is a classical triorthogonal space (NOTES Def 1, Thm 1).

    ``basis`` is a 2D array of GF(p) rows b_1..b_K.  Checks the quadratic condition
    sum_i (b_a*b_b)_i == 0 for every a <= b AND the cubic condition
    sum_i (b_a*b_b*b_c)_i == 0 for every a <= b <= c, all mod p, **including repeated
    indices** (a==b, a==b==c, ...).  If the basis rows satisfy these for all unordered
    index multisets, every triple of vectors in the span does too (the star product is
    multilinear), so the space is triorthogonal.
    """
    assert isinstance(p, int) and p >= 2, "p must be an integer prime >= 2"
    B = np.asarray(basis, dtype=object) % p  # python ints, exact mod-p arithmetic
    B = np.array(B.tolist(), dtype=object)
    K = B.shape[0]

    # Quadratic (pair) condition, including a == b.
    for a in range(K):
        for b in range(a, K):
            s = int(np.sum(B[a] * B[b])) % p
            if s != 0:
                return False

    # Cubic (triple) condition, including repeated indices.
    for a in range(K):
        for b in range(a, K):
            ab = B[a] * B[b]
            for c in range(b, K):
                s = int(np.sum(ab * B[c])) % p
                if s != 0:
                    return False
    return True


def puncture_matrix(G, columns0, p):
    """Delete the given 0-indexed columns from G -> punctured generator Gp (NOTES sec 3, eq 1).

    Puncturing removes the chosen coordinates; the resulting generator G' spans the
    punctured code PRM (NOTES Thm 3).  Returns a GF(p) array with ``len(columns0)``
    fewer columns.
    """
    GF = GFp(p)
    G = GF(np.asarray(G) % p)
    cols0 = list(columns0)
    keep = np.ones(G.shape[1], dtype=bool)
    if cols0:
        keep[np.asarray(cols0, dtype=np.int64)] = False
    return G[:, keep]


def shorten_matrix(G, columns0, p):
    """Generator of the shortened code G0 (NOTES sec 3, eq 1; Thm 3 duality).

    The shortened code keeps only codewords of rowspan(G) that are zero on ``columns0``,
    then deletes those columns.  A combination u @ G is zero on columns0 iff
    u @ G[:, columns0] == 0, i.e. u is in the left null space of the punctured-column
    submatrix; those combinations, restricted to the surviving columns, give G0.  By
    the Thm 3 duality SRM(r,m;S)^perp = PRM(rtilde,m;S), so dual(G0) = G'^perp at the
    dual degree.
    """
    GF = GFp(p)
    G = GF(np.asarray(G) % p)
    cols0 = list(columns0)
    if not cols0:
        return G
    sub = G[:, np.asarray(cols0, dtype=np.int64)]  # (K, k): the punctured columns
    u = sub.left_null_space()                       # rows u with u @ sub == 0
    if u.shape[0] == 0:
        return GF.Zeros((0, G.shape[1] - len(cols0)))
    shortened = u @ G                               # codewords zero on columns0
    return puncture_matrix(shortened, cols0, p)


def dual_matrix(G, p):
    """Generator of the dual (null space) code of G via galois (NOTES sec 3; Lemma 1).

    The dual code is { y : G y^T = 0 } whose basis is the (right) null space of G.
    Used for the Z-stabilizers G'^perp of the CSS code.
    """
    GF = GFp(p)
    G = GF(np.asarray(G) % p)
    return G.null_space()


def build_triorthogonal_code(p, m, r, puncture_columns_1indexed, G=None) -> dict:
    """MASTER builder: the CSS(G0 -> X, G'^perp -> Z) triorthogonal code (NOTES sec 3, sec 10 item 7).

    Builds G = rm_generator(r,m,p) (a triorthogonal space, requiring 3r < m(p-1), Thm 1),
    punctures and shortens at the given 1-indexed columns (mapped to 0-indexed via c-1,
    NOTES App. C), and assembles the CSS code:

        X_stab = G0                       (shortened generator)
        Z_stab = dual(Gp) = G'^perp       (NOTES Thm 3: = PRM(rtilde,m;S))

    Parameters of the [[n, k, d]]_p code: n = p**m - #punctures, k = #punctures.  This
    labeling is valid only when the surviving (punctured-column) submatrix of G is full
    rank == #punctures (NOTES sec 3 caveat, sec 10 item 9); ``full_rank`` reports this.

    Returns dict with keys: n, k, full_rank (bool), Gp, G0, Z_stab, X_stab.
    """
    from .reedmuller import rm_generator

    assert isinstance(p, int) and p >= 2, "p must be an integer prime >= 2"
    assert isinstance(m, int) and m >= 1, "m must be a positive integer"
    assert 3 * r < m * (p - 1), "RM space is triorthogonal only when 3r < m(p-1) (Thm 1)"

    GF = GFp(p)
    if G is None:  # callers (e.g. the search) may pass a prebuilt generator to reuse it
        G = rm_generator(r, m, p)

    cols0 = [c - 1 for c in puncture_columns_1indexed]
    npunc = len(cols0)
    n0 = p ** m
    assert all(0 <= c < n0 for c in cols0), "puncture column out of range [1, p**m]"
    assert len(set(cols0)) == npunc, "duplicate puncture columns"

    Gp = puncture_matrix(G, cols0, p)
    G0 = shorten_matrix(G, cols0, p)

    # Full-rank check: the punctured-column submatrix must have rank == #punctures,
    # otherwise the [[n,k,d]] labeling is invalid (NOTES sec 3 caveat).
    if npunc == 0:
        full_rank = True
    else:
        sub = GF(np.asarray(G[:, np.asarray(cols0, dtype=np.int64)]) % p)
        full_rank = int(np.linalg.matrix_rank(sub)) == npunc

    Z_stab = dual_matrix(Gp, p)
    X_stab = G0

    return {
        "n": n0 - npunc,
        "k": npunc,
        "full_rank": bool(full_rank),
        "Gp": Gp,
        "G0": G0,
        "Z_stab": Z_stab,
        "X_stab": X_stab,
    }
