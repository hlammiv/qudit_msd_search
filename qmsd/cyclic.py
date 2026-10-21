"""Cyclic triorthogonal codes from spectral supports (arXiv:2610.08012 secs 2, 4).

A cyclic code C of length n over F_p (p does not divide n) is fixed by its spectral support
S, a union of p-cyclotomic cosets mod n: C = <g> with g(x) = prod_{j not in S} (x - alpha^j)
for alpha a primitive n-th root of unity in GF(p^m), m = ord_n(p).  Any choice of alpha
changes S by a multiplier, giving a permutation-equivalent code, so weight distributions
(and hence distillation performance) do not depend on it.

Two constructions (paper Props 12, 13):

  * A: stabilizers = basis of C, logical row = all-ones.  Admissible iff 0 not in S, 2S, 3S.
  * B: W = C is a triorthogonal space iff 0 not in 3S; puncture W at coordinate 0 (Lemma 11)
       to get an [[n-1, 1, d(W^perp)-1]] code with lambda = -1.

This module is a minimal builder used as a validation oracle for ``qmsd.performance``; it
returns ``(G0, H1)`` = (stabilizer rows, logical rows) as int64 arrays mod p.
"""
from __future__ import annotations

from functools import lru_cache

import galois
import numpy as np


def multiplicative_order(p: int, n: int) -> int:
    """Least m >= 1 with n | p^m - 1 (requires gcd(n, p) = 1)."""
    assert n >= 1 and n % p != 0, "length must be coprime to p"
    m, x = 1, p % n
    while x != 1 % n:
        x = (x * p) % n
        m += 1
    return m


def cyclotomic_cosets(p: int, n: int) -> list[tuple[int, ...]]:
    """The p-cyclotomic cosets {s, ps, p^2 s, ...} mod n, each sorted, ordered by min element."""
    seen, out = set(), []
    for s in range(n):
        if s in seen:
            continue
        coset, x = [], s
        while x not in coset:
            coset.append(x)
            x = (x * p) % n
        seen.update(coset)
        out.append(tuple(sorted(coset)))
    return out


def sumset(S, t: int, n: int) -> set[int]:
    """tS = S + ... + S (t summands) mod n."""
    out = {0}
    for _ in range(t):
        out = {(a + s) % n for a in out for s in S}
    return out


def admissible_A(S, n: int) -> bool:
    """Construction A criterion: 0 not in S, 2S, 3S (Prop 12)."""
    return all(0 not in sumset(S, t, n) for t in (1, 2, 3))


def admissible_B(S, n: int) -> bool:
    """Construction B criterion: 0 not in 3S (Prop 13)."""
    return 0 not in sumset(S, 3, n)


@lru_cache(maxsize=None)
def _root_field(p: int, n: int):
    m = multiplicative_order(p, n)
    GF = galois.GF(p ** m)
    beta = GF.primitive_element ** ((p ** m - 1) // n)
    return GF, beta


def cyclic_generator(p: int, n: int, S) -> np.ndarray:
    """Generator-polynomial coefficients g_0..g_{n-|S|} (ascending, in F_p) of the code with support S."""
    S = {s % n for s in S}
    assert S == {(p * s) % n for s in S}, "S must be a union of p-cyclotomic cosets"
    GF, beta = _root_field(p, n)
    g = galois.Poly([1], field=GF)
    for j in range(n):
        if j not in S:
            g = g * galois.Poly([1, -(beta ** j)], field=GF)
    coeffs = [int(c) for c in g.coeffs[::-1]]  # ascending powers
    assert all(0 <= c < p for c in coeffs), "generator not over the prime field (S not coset-closed?)"
    return np.array(coeffs, dtype=np.int64)


def cyclic_basis(p: int, n: int, S) -> np.ndarray:
    """|S| x n basis of C = <g>: the shifts x^i g(x), i < |S|."""
    g = cyclic_generator(p, n, S)
    dim = n - (len(g) - 1)
    rows = np.zeros((dim, n), dtype=np.int64)
    for i in range(dim):
        rows[i, i:i + len(g)] = g
    return rows % p


def construction_A(p: int, n: int, S) -> tuple[np.ndarray, np.ndarray]:
    """(G0, H1) for construction A: stabilizers = basis of C, logical row = all-ones (Prop 12)."""
    assert admissible_A(S, n), "support not admissible for construction A"
    return cyclic_basis(p, n, S), np.ones((1, n), dtype=np.int64)


def construction_B(p: int, n: int, S) -> tuple[np.ndarray, np.ndarray]:
    """(G0, H1) for construction B: W = C punctured at coordinate 0 (Prop 13, Lemma 11).

    Length of the returned code is n - 1.  The logical row is a codeword g with g_0 = 1;
    the stabilizers are a basis of {w in W : w_0 = 0}, both with coordinate 0 deleted.
    """
    assert admissible_B(S, n), "support not admissible for construction B"
    W = cyclic_basis(p, n, S)
    piv = int(np.flatnonzero(W[:, 0])[0])
    g = W[piv] * pow(int(W[piv, 0]), -1, p) % p
    others = [(W[i] - W[i, 0] * g) % p for i in range(W.shape[0]) if i != piv]
    G0 = np.array(others, dtype=np.int64).reshape(-1, n)[:, 1:]
    return G0 % p, g[None, 1:] % p
