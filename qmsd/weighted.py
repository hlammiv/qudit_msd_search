"""Weighted Reed-Muller supports: shorter triorthogonal spaces from charge distributions.

In the charge picture of arXiv:2610.29457 sec 4, a set of points X in F_p^m with charges mu supports
RM(r, m) as a triorthogonal space (with coefficients mu) iff the charges annihilate every polynomial
of degree <= 3r:

    sum_{x in X} mu_x f(x) = 0    for deg f <= 3r.

The unit charges on all of F_p^m are the usual punctured-RM setting (3r < m(p-1)).  More generally,
mu can be the values of ANY polynomial g of degree <= m(p-1) - 1 - 3r, i.e. a codeword of
RM(3r, m)^perp, and then X = supp(g) can be much smaller than p^m.  Example: p = 7, m = 3, r = 5
leaves degree 2, and g = x1 (x1 - 1) gives |X| = 245 instead of 343.

Puncturing k points of X gives a k-output code as usual.  Output i implements U^{mu_i}, so outputs
should sit at points whose charge is a cube (every unit for p = 3, 5 mod 6; for p = 1 mod 3 only a
third of them).  Inputs at non-cube charges consume noisy U^{mu} magic states, which arXiv:2610.29457
counts as one input each (its [[20,1,4]]_7 does the same).
"""
from __future__ import annotations

import numpy as np

from .reedmuller import points, r_max, rm_generator
from .triorthogonal import is_cube_mod_p


def charge_degree_budget(p, m, r) -> int:
    """Largest degree of a charge polynomial g that keeps RM(r, m) triorthogonal on supp(g)."""
    return m * (p - 1) - 1 - 3 * r


def weighted_support(p, m, g_values, r=None):
    """(G restricted to X, X, mu) for charges mu = g(x) on X = supp(g).  ``g_values``: g at the
    p^m points in ``reedmuller.points`` order (the caller guarantees deg g <= the budget)."""
    r = r_max(m, p) if r is None else r
    assert charge_degree_budget(p, m, r) >= 0, "no slack: only the full space with unit charges works"
    g = np.asarray(g_values, dtype=np.int64) % p
    X = np.nonzero(g)[0]
    G = np.asarray(rm_generator(r, m, p), dtype=np.int64) % p
    return G[:, X], X, g[X]


def weighted_code(GX, mu, puncture, p):
    """(G0, H1, coeffs) after puncturing columns ``puncture`` (indices into X) of the weighted space.

    H1 rows are the canonical logical rows (delta on the punctures), G0 the shortened space; the
    surviving coefficients are ``coeffs``.  Asserts the punctured columns are independent.
    """
    from .performance import canonical_logical_rows
    from .triorthogonal import shorten_matrix

    puncture = list(puncture)
    keep = [l for l in range(GX.shape[1]) if l not in set(puncture)]
    G0 = np.asarray(shorten_matrix(GX, puncture, p), dtype=np.int64) % p
    H1 = canonical_logical_rows(GX, puncture, p)
    H1 = np.asarray(H1, dtype=np.int64) % p
    return G0, H1, np.asarray(mu, dtype=np.int64)[keep] % p


def cube_points(mu, p):
    """Indices (into X) whose charge is a cube, i.e. where an output stays in U's Clifford class."""
    return [l for l, u in enumerate(np.asarray(mu) % p) if is_cube_mod_p(int(u), p)]


def quadric(p, m, diag, lin=(), const=0):
    """g(x) = sum_i diag[i] x_i^2 + sum_i lin[i] x_i + const at all points.  Up to affine
    equivalence every degree-2 charge is one of these (diagonalise the quadratic part).

    At p = 7, m = 3 the best-spread output classes come from the parabolic quadric
    ``quadric(7, 3, (1, 1, 0), (0, 0, 1))`` and from ``quadric(7, 3, (1, 1, 1), const=1)``:
    at most 16 class points on a plane and 4 on a line, against 49 on a plane for x1(x1-1).
    """
    P = np.asarray(points(m, p), dtype=np.int64)
    diag = np.asarray(diag, dtype=np.int64)
    lin = np.asarray(list(lin) + [0] * (m - len(lin)), dtype=np.int64)
    return (P * P @ diag + P @ lin + const) % p


def class_points(mu, p, unit):
    """Indices (into X) whose charge is ``unit`` times a cube; ``unit=1`` is ``cube_points``.
    Multiplying g by a non-cube unit relabels the classes, so outputs may be placed in any one
    class and the whole code rescaled afterwards."""
    inv = pow(int(unit), -1, p)
    return [l for l, u in enumerate(np.asarray(mu) % p) if is_cube_mod_p(int(u) * inv % p, p)]


def product_of_linear(p, m, factors):
    """g(x) = prod (a . x - b) over ``factors`` [(a, b), ...], evaluated at all points."""
    P = np.asarray(points(m, p), dtype=np.int64)
    g = np.ones(len(P), dtype=np.int64)
    for a, b in factors:
        g = g * ((P @ np.asarray(a, dtype=np.int64) - b) % p) % p
    return g
