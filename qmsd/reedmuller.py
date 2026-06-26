"""Generalized Reed-Muller codes over F_p (NOTES sec 4; Thm 1, Thm 2).

A Reed-Muller code RM_p(r,m) is the evaluation code (NOTES eq 2)

    RM_p(r,m) = { (f(v) : v in F_p^m) : deg f <= r,  f in F_p[x_1..x_m]/<x_i^p - x_i> }.

Because x_i^p = x_i in F_p, each variable's individual degree is at most p-1, so the
monomial basis is { x_1^{a_1}...x_m^{a_m} : 0 <= a_i <= p-1, sum a_i <= r }.  Each such
monomial, evaluated over all p^m points, is one generator row.

Global conventions (NOTES sec 4 "Typo warning"): the triorthogonality / dual degrees
use m(p-1), NOT the paper's misprinted p(m-1).
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product

import numpy as np

from .field import GFp


def r_max(m, p) -> int:
    """floor((m(p-1)-1)/3): largest r with 3r < m(p-1) (NOTES Thm 1, eq near 127).

    Uses the corrected m(p-1) threshold (NOTES sec 4 typo warning), matching the
    oracle's OracleCode.r_max.
    """
    assert isinstance(p, int) and p >= 2, "p must be an integer prime >= 2"
    assert isinstance(m, int) and m >= 1, "m must be a positive integer"
    r = (m * (p - 1) - 1) // 3
    assert 3 * r < m * (p - 1), "r_max must satisfy 3r < m(p-1)"
    return r


def r_tilde(r, m, p) -> int:
    """m(p-1) - r - 1: the dual degree (NOTES Lemma 1, eq 113).

    RM_p(r,m)^perp = RM_p(rtilde, m); rtilde controls the quantum distance.
    """
    assert isinstance(p, int) and p >= 2, "p must be an integer prime >= 2"
    assert isinstance(m, int) and m >= 1, "m must be a positive integer"
    assert 0 <= r <= m * (p - 1), "r must lie in [0, m(p-1)]"
    return m * (p - 1) - r - 1


def d_rm(r, m, p) -> int:
    """Minimum distance of RM_p(r,m) via Schwartz-Zippel (NOTES Thm 2, eqs 4-5).

    d_RM(r,m) = p^(m - floor(r/(p-1))) * (1 - (r mod (p-1))/p), computed exactly in
    integers as  p^(m - r//(p-1)) * (p - r%(p-1)) // p.
    """
    assert isinstance(p, int) and p >= 2, "p must be an integer prime >= 2"
    assert isinstance(m, int) and m >= 1, "m must be a positive integer"
    assert 0 <= r <= m * (p - 1), "r must lie in [0, m(p-1)]"
    alpha = r // (p - 1)
    beta = r % (p - 1)
    return p ** (m - alpha) * (p - beta) // p


def monomial_exponents(m, r, p) -> list:
    """All exponent tuples (a_1..a_m), 0<=a_i<=p-1, sum(a_i)<=r (NOTES sec 4, line 444).

    Deterministic order: lexicographic over itertools.product(range(p), repeat=m)
    (a_1 varying fastest), filtered to total degree <= r.
    """
    assert isinstance(p, int) and p >= 2, "p must be an integer prime >= 2"
    assert isinstance(m, int) and m >= 1, "m must be a positive integer"
    assert 0 <= r <= m * (p - 1), "r must lie in [0, m(p-1)]"
    return [a for a in product(range(p), repeat=m) if sum(a) <= r]


@lru_cache(maxsize=None)
def points(m, p):
    """Evaluation points: ndarray (p**m, m); row c = base-p digits of c (NOTES App. C).

    x_1 is least significant: point c = (c % p, (c // p) % p, ..., (c // p^{m-1}) % p),
    for c = 0 .. p**m - 1.  (Column convention is 1-indexed; column c <-> point c-1.)
    """
    assert isinstance(p, int) and p >= 2, "p must be an integer prime >= 2"
    assert isinstance(m, int) and m >= 1, "m must be a positive integer"
    n = p ** m
    pts = np.empty((n, m), dtype=np.int64)
    for j in range(m):
        pts[:, j] = (np.arange(n) // (p ** j)) % p
    return pts


def rm_generator(r, m, p):
    """Generator of RM_p(r,m): GF(p) array (n_monomials, p**m) (NOTES sec 4, eq 2).

    Row i is monomial i = prod_j x_j^{a_j} (the i-th tuple from monomial_exponents)
    evaluated over every point of points(m,p).  The 0^0 = 1 convention is honored by
    treating any a_j == 0 factor as the constant 1 (so it contributes even at x_j = 0).
    """
    GF = GFp(p)
    exps = monomial_exponents(m, r, p)
    pts = GF(points(m, p).astype(np.int64) % p)  # (p**m, m) over F_p
    n = p ** m
    G = GF.Zeros((len(exps), n))
    for i, a in enumerate(exps):
        row = GF.Ones(n)
        for j, aj in enumerate(a):
            if aj:
                row = row * (pts[:, j] ** aj)
        G[i, :] = row
    return G
