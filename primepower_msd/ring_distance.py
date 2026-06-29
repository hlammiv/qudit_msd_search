"""Exact (bounded-weight) minimum distance for ring CSS codes — the M3 distance engine.

distance = min Hamming weight of a nontrivial logical operator:
  * X-logical : v in M_Z^perp \\ M_X   (v . b == 0 for all Z-stabilizers b;  v not in M_X)
  * Z-logical : v in M_X^perp \\ M_Z   (v . a == 0 for all X-stabilizers a;  v not in M_Z)

We enumerate weight-1, weight-2, ... vectors (cheap for small weight) and return the first weight at
which a nontrivial logical appears.  Scalable: uses only generator-level membership (the scalable
Howell kernel / in_span), never the d^n grid.  Returns `cap + 1` (a lower bound) if distance > cap.
"""

from __future__ import annotations

from itertools import combinations, product

from .ringlinalg import in_span


def _ortho(v, gens, d) -> bool:
    """v . b == 0 (mod d) for every generator row b (i.e. v in <gens>^perp)."""
    return all(sum(v[i] * int(b[i]) for i in range(len(v))) % d == 0 for b in gens)


def min_logical_weight(mx_gens, mz_gens, n: int, d: int, cap: int = 6) -> int:
    """Minimum weight of a nontrivial X- or Z-logical operator, searched up to `cap` (else cap+1)."""
    mxg = [list(map(int, r)) for r in mx_gens]
    mzg = [list(map(int, r)) for r in mz_gens]
    for w in range(1, cap + 1):
        for pos in combinations(range(n), w):
            for vals in product(range(1, d), repeat=w):
                v = [0] * n
                for p, val in zip(pos, vals):
                    v[p] = val
                # X-logical: in M_Z^perp, not in M_X
                if _ortho(v, mzg, d) and not in_span(v, mxg, d):
                    return w
                # Z-logical: in M_X^perp, not in M_Z
                if _ortho(v, mxg, d) and not in_span(v, mzg, d):
                    return w
    return cap + 1


def code_gamma(n: int, k: int, distance: int) -> float:
    """Overhead exponent gamma = log(n/k)/log(distance).  gamma < 1 == sublogarithmic (needs k>=2)."""
    import math
    if distance <= 1 or k < 1:
        return float("inf")
    return math.log(n / k) / math.log(distance)
