"""Qudit doubling (Prakash & Ramya, arXiv:2610.29457): short single-output triorthogonal codes.

Theorem 4 doubles a triorthogonal seed T = (y; T0) of length n_T with a (weighted) P-code S of
length n_P -- a self-orthogonal code orthogonal to 1 with total charge N = sum(mu) != 0 -- into a

    [[2 n_P + n_T, 1, min(d_P, d_T + 2)]]_p

code with the seed's logical exponent.  Rows (columns: copy 1 | copy 2 | seed):

    L = (-k1,  k1, y)          logical
        ( s,   -s, 0)          s in a basis of S
    K = ( 0,   k1, y)          linking row
        ( 0,    0, t)          t in a basis of T0

with coefficients eta*mu on both copies and the seed's own on the seed, where eta*N*k^3 + lambda_T = 0.

The P-codes of the paper's Table 1 that its doubled codes use are built here, and
``paper_codes(p)`` returns the doubled codes of its Table 3 as (label, G0, H1, coeffs).  Under
twirled noise the coefficients do not enter the error model, only the lambda (Clifford class) check.
"""
from __future__ import annotations

import itertools
from functools import lru_cache

import numpy as np

from .triorthogonal import logical_lambdas


def _rows(*strs):
    return np.array([[int(ch) for ch in s] for s in strs], dtype=np.int64)


def double(S, seed_H, p, mu=None, seed_c=None, eta=None, kappa=None):
    """Theorem 4.  ``S``: r x n_P basis of the P-code (charges ``mu``, default 1).  ``seed_H``: the
    seed's rows, logical row first (coefficients ``seed_c``, default 1).  Returns (H, coeffs) with
    the logical row first.  ``eta``/``kappa`` default to the first solution of eta*N*kappa^3 + lambda_T = 0.
    """
    S = np.atleast_2d(np.asarray(S, dtype=np.int64) % p)
    T = np.atleast_2d(np.asarray(seed_H, dtype=np.int64) % p)
    nP, nT = S.shape[1], T.shape[1]
    mu = np.ones(nP, dtype=np.int64) if mu is None else np.asarray(mu, dtype=np.int64) % p
    seed_c = np.ones(nT, dtype=np.int64) if seed_c is None else np.asarray(seed_c, dtype=np.int64) % p
    N = int(mu.sum() % p)
    lam = logical_lambdas(T, 1, p, seed_c)[0]
    assert N != 0 and lam != 0, "P-code needs total charge != 0 and the seed lambda != 0"
    if eta is None or kappa is None:
        eta, kappa = next((e, k) for e in range(1, p) for k in range(1, p)
                          if (e * N * k ** 3 + lam) % p == 0)
    assert (eta * N * kappa ** 3 + lam) % p == 0, "eta * N * kappa^3 + lambda_T must vanish"
    one, y = np.ones(nP, dtype=np.int64), T[0]
    rows = [np.concatenate([-kappa * one, kappa * one, y])]
    rows += [np.concatenate([s, -s, np.zeros(nT, dtype=np.int64)]) for s in S]
    rows.append(np.concatenate([0 * one, kappa * one, y]))
    rows += [np.concatenate([0 * one, 0 * one, t]) for t in T[1:]]
    coeffs = np.concatenate([eta * mu, eta * mu, seed_c]) % p
    return np.array(rows, dtype=np.int64) % p, coeffs


def is_pcode(S, p, mu=None) -> bool:
    """Definition 2 / 10: sum mu s s' = 0 and sum mu s = 0 on S, and N = sum mu != 0."""
    S = np.atleast_2d(np.asarray(S, dtype=np.int64) % p)
    mu = np.ones(S.shape[1], dtype=np.int64) if mu is None else np.asarray(mu, dtype=np.int64) % p
    return (int(mu.sum() % p) != 0 and not np.any((S * mu) @ S.T % p)
            and not np.any((S * mu).sum(axis=1) % p))


def pcode_distance(S, p, max_words=5_000_000) -> int:
    """d_P = min wt(v), v in S^perp (plain dot product) with sum(v) != 0, by enumerating S^perp."""
    import galois

    GF = galois.GF(p)
    Sp = np.asarray(GF(np.atleast_2d(np.asarray(S) % p)).null_space(), dtype=np.int64)
    assert p ** Sp.shape[0] <= max_words, "S^perp too large to enumerate"
    best = None
    for msg in itertools.product(range(p), repeat=Sp.shape[0]):
        v = np.asarray(msg, dtype=np.int64) @ Sp % p
        if v.sum() % p:
            w = int(np.count_nonzero(v))
            best = w if best is None else min(best, w)
    return best


@lru_cache(maxsize=None)
def cyclic_pcode(p: int, n: int, sigma: frozenset, d: int) -> np.ndarray:
    """Basis of the cyclic P-code with spectral support ``sigma`` (Lemma 15) and distance ``d``.

    The paper fixes its n-th root of unity via SymPy; ours differs by a multiplier t, so try the
    supports t*sigma (t coprime to n) until the P-code conditions hold with distance ``d``.
    """
    from math import gcd

    from .cyclic import cyclic_basis

    for t in [t for t in range(1, n) if gcd(t, n) == 1]:
        S = cyclic_basis(p, n, {(t * s) % n for s in sigma})
        if is_pcode(S, p) and pcode_distance(S, p) == d:
            return S
    raise ValueError(f"no multiplier of the support gives a [[{n},1,{d}]]_{p} P-code")


def quadratic_residues(q: int) -> frozenset:
    return frozenset((x * x) % q for x in range(1, q))


# --- the P-codes of Table 1 used by the doubled codes ---------------------------------------
PCODE_7W_3 = (np.array([(0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2)]).T,
              np.array([2, 2, 1, 2, 1, 1, 2]))                       # Example 11, [[7,1,3]]^P_3
PCODE_10_3 = _rows("1000210122", "0100221012", "0010222101", "0001212210")   # [[10,1,4]]^P_3
PCODE_7W_7 = (_rows("0021114", "0401236", "2001323"), np.array([2, 2, 2, 4, 4, 1, 2]))  # Ex. 14


def paper_codes(p: int) -> list[tuple[str, np.ndarray, np.ndarray, np.ndarray]]:
    """The doubled codes of arXiv:2610.29457 Table 3 that are cheap to enumerate exactly:
    (label, G0, H1, coeffs).  [[72,1,6]]_3 (18 stabilizers) and [[118,1,8]]_3 are left out."""
    from .cyclic import construction_B
    from .protocol import _rm_k1

    out = []

    def add(label, H, c):
        out.append((label, H[1:], H[:1], c))

    if p == 3:
        trivial = np.array([[1]])
        H15, c15 = double(*PCODE_7W_3[:1], trivial, 3, mu=PCODE_7W_3[1])
        add("[[15,1,3]]_3 dbl", H15, c15)
        G0, H1 = _rm_k1(3, 2)                                            # [[8,1,2]]_3 seed
        H28, c28 = double(PCODE_10_3, np.vstack([H1, G0]), 3)
        add("[[28,1,4]]_3 dbl", H28, c28)
        golay = cyclic_pcode(3, 11, quadratic_residues(11), 5)          # [[11,1,5]]^P_3
        H37, c37 = double(golay, H15, 3, seed_c=c15)
        add("[[37,1,5]]_3 dbl", H37, c37)
    elif p == 5:
        G0, H1 = construction_B(5, 12, {2, 3, 9, 10})                   # cyclic [[11,1,3]]_5 seed
        qr = cyclic_pcode(5, 11, quadratic_residues(11), 5)             # [[11,1,5]]^P_5
        H33, c33 = double(qr, np.vstack([H1, G0]), 5)
        add("[[33,1,5]]_5 dbl", H33, c33)
    elif p == 7:
        G0, H1 = _rm_k1(7, 1)                                            # [[6,1,2]]_7 RS seed
        H20, c20 = double(PCODE_7W_7[0], np.vstack([H1, G0]), 7, mu=PCODE_7W_7[1], eta=2, kappa=3)
        add("[[20,1,4]]_7 dbl", H20, c20)
    return out


def prakash_saha_d2(m: int) -> tuple[np.ndarray, np.ndarray]:
    """(G0, H1) of the [[6m+2, 3m-2, 2]]_3 code of arXiv:2403.06228: the 9m-point triorthogonal
    space spanned by w = (0,1,2)^{3m} and v^(a) (1s on block a, 2s on the last block), punctured
    at the first trit of blocks 0..3m-3."""
    n, k = 9 * m, 3 * m - 2
    T = [np.tile([0, 1, 2], 3 * m)]
    for a in range(1, 3 * m):
        v = np.zeros(n, dtype=np.int64)
        v[3 * (a - 1):3 * a] = 1
        v[n - 3:] = 2
        T.append(v)
    T = np.array(T, dtype=np.int64)
    keep = [c for c in range(n) if c not in {3 * j for j in range(k)}]
    return np.array([T[0], T[k + 1]])[:, keep] % 3, T[1:k + 1][:, keep] % 3
