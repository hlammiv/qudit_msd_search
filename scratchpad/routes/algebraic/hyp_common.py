"""Shared helpers: evaluate a random polynomial of degree <= delta on F_p^m, zero set,
rank of RM(r,m) on a point set, and the theorem bound d(RM(rtilde+delta, m))."""
import numpy as np
from qmsd.reedmuller import rm_generator, points, monomial_exponents, d_rm
from qmsd.field import GFp


def eval_poly(coeffs, delta, m, p):
    """coeffs: array in monomial_exponents(m, delta, p) order -> values on points(m, p)."""
    G = np.asarray(rm_generator(delta, m, p), dtype=np.int64)
    c = np.asarray(coeffs, dtype=np.int64) % p
    return (c @ G) % p


def random_poly_values(delta, m, p, rng, homogeneous_top=False):
    exps = monomial_exponents(m, delta, p)
    c = rng.integers(0, p, size=len(exps))
    if homogeneous_top:  # force a nonzero top-degree part
        top = [i for i, a in enumerate(exps) if sum(a) == delta]
        if (c[top] == 0).all():
            c[top[0]] = 1
    return eval_poly(c, delta, m, p), c


def rank_modp(M, p):
    GF = GFp(p)
    return int(np.linalg.matrix_rank(GF(np.asarray(M, dtype=np.int64) % p)))


def bound(delta, r, m, p):
    """Lower bound on the weight of kernel words meeting F^m \\ Z(g), deg g = delta:
    d(RM(rtilde + delta, m)) with rtilde = m(p-1) - r - 1."""
    rt = m * (p - 1) - r - 1
    return d_rm(rt + delta, m, p) if rt + delta <= m * (p - 1) else None
