"""Galois ring GR(4, m) = Z_4[x]/(h(x)), h a Hensel-lifted primitive polynomial, and the non-CSS
code builders that the Phase-5 gamma<1 codes came from (generalizes the untracked gr_engine.py, which
was hard-wired to GR(4,2) = Z_4[xi]/(xi^2 + xi + 1)).

Elements are length-m int tuples (coefficients of 1, x, ..., x^{m-1}) mod 4.  Additive group (Z_4)^m,
so one GR coordinate = m genuine order-4 ring qudits (NOT the GF(4) collapse).

  * `hensel_lift_primitive(h2)`: Graeffe lift of a primitive binary polynomial to Z_4 so that the lifted
    root xi has order 2^m - 1 (Teichmueller set T = {0, 1, xi, ..., xi^{2^m-2}}).
  * Frobenius sigma(f(xi)) = f(xi^2) (ring automorphism; trace Tr(u) = sum_j sigma^j(u) in Z_4).
  * Trace Gram G_ij = Tr(b_i b_j) on the monomial basis is NON-DEGENERATE (det a unit) -- the lever that
    makes the GR Reed-Muller twist a genuine non-CSS code (CONSTRUCTIVE_PLAN.md).
  * Physical embedding: X-part of a GR coordinate u = its basis components; Z-part of v = G @ components(v),
    so the STANDARD symplectic form equals the trace-symplectic form Tr(sum u_i v'_i - u'_i v_i).

For m = 2 this reproduces gr_engine exactly (h = x^2 + x + 1, GRAM = [[2,3],[3,3]]).
"""

from __future__ import annotations

import itertools as it
from functools import lru_cache

import numpy as np

# primitive binary polynomials (coefficient lists low -> high), one per degree
PRIMITIVE_F2 = {
    1: [1, 1],                 # x + 1
    2: [1, 1, 1],              # x^2 + x + 1
    3: [1, 1, 0, 1],           # x^3 + x + 1
    4: [1, 1, 0, 0, 1],        # x^4 + x + 1
    5: [1, 0, 1, 0, 0, 1],     # x^5 + x^2 + 1
    6: [1, 1, 0, 0, 0, 0, 1],  # x^6 + x + 1
}


def _polymul(a, b, mod):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % mod
    return out


def hensel_lift_primitive(h2):
    """Graeffe lift: h4(x^2) = +-(e(x)^2 - o(x)^2) mod 4, with h2 = e + o (even/odd parts).
    The result is a monic degree-m polynomial over Z_4 reducing to h2 mod 2 whose root has order 2^m - 1."""
    m = len(h2) - 1
    e = [c if i % 2 == 0 else 0 for i, c in enumerate(h2)]
    o = [c if i % 2 == 1 else 0 for i, c in enumerate(h2)]
    e2 = _polymul(e, e, 4)
    o2 = _polymul(o, o, 4)
    diff = [(x - y) % 4 for x, y in zip(e2 + [0] * (len(o2) - len(e2)), o2 + [0] * (len(e2) - len(o2)))]
    # diff is a polynomial in x^2: take even coefficients
    h4 = [diff[2 * i] for i in range(m + 1)]
    if h4[m] % 4 == 3:        # make monic
        h4 = [(-c) % 4 for c in h4]
    assert h4[m] % 4 == 1
    return h4


class GR:
    """GR(4, m) with element tuples of length m (mod 4)."""

    def __init__(self, m: int):
        self.m = m
        self.h = hensel_lift_primitive(PRIMITIVE_F2[m])          # monic, degree m
        # reduction of x^m: x^m = -sum_{i<m} h_i x^i
        self.xm = [(-c) % 4 for c in self.h[:m]]
        self.zero = tuple([0] * m)
        self.one = tuple([1] + [0] * (m - 1))
        self.xi = tuple([0, 1] + [0] * (m - 2)) if m >= 2 else tuple([3])   # for m=1: x = -1 = 3
        self.elems = [tuple(t) for t in it.product(range(4), repeat=m)]
        self.basis = [tuple(1 if j == i else 0 for j in range(m)) for i in range(m)]
        self.gram = np.array([[self.trace(self.mul(bi, bj)) for bj in self.basis] for bi in self.basis], dtype=np.int64) % 4
        self.gram_det_unit = int(round(np.linalg.det(self.gram))) % 2 == 1
        self._teich = None

    # ---- arithmetic
    def add(self, u, v):
        return tuple((a + b) % 4 for a, b in zip(u, v))

    def neg(self, u):
        return tuple((-a) % 4 for a in u)

    def smul(self, c: int, u):
        return tuple((c * a) % 4 for a in u)

    def mul(self, u, v):
        m = self.m
        prod = [0] * (2 * m - 1)
        for i, a in enumerate(u):
            if a:
                for j, b in enumerate(v):
                    prod[i + j] = (prod[i + j] + a * b) % 4
        # reduce degrees >= m
        for deg in range(2 * m - 2, m - 1, -1):
            c = prod[deg]
            if c:
                prod[deg] = 0
                for i in range(m):
                    prod[deg - m + i] = (prod[deg - m + i] + c * self.xm[i]) % 4
        return tuple(prod[:m])

    def pow(self, u, e: int):
        r = self.one
        base = u
        while e:
            if e & 1:
                r = self.mul(r, base)
            base = self.mul(base, base)
            e >>= 1
        return r

    def frobenius(self, u):
        """sigma(f(x)) = f(x^2): ring automorphism of GR(4,m) (x is the Teichmueller root)."""
        out = self.zero
        for i, a in enumerate(u):
            if a:
                out = self.add(out, self.smul(a, self.pow(self.xi, 2 * i)))
        return out

    def trace(self, u) -> int:
        s = self.zero
        v = u
        for _ in range(self.m):
            s = self.add(s, v)
            v = self.frobenius(v)
        assert all(c == 0 for c in s[1:]), "trace must land in Z_4"
        return int(s[0]) % 4

    def is_unit(self, u) -> bool:
        return any(c % 2 for c in u)

    def teichmueller(self):
        """T = {0} U {xi^j : j < 2^m - 1}: the multiplicatively closed lift of GF(2^m)."""
        if self._teich is None:
            T = [self.zero]
            p = self.one
            for _ in range(2 ** self.m - 1):
                T.append(p)
                p = self.mul(p, self.xi)
            assert p == self.one, "xi must have order 2^m - 1"
            self._teich = T
        return self._teich

    # ---- physical embedding (standard symplectic == trace-symplectic)
    def xpart(self, u) -> np.ndarray:
        return np.array(u, dtype=np.int64) % 4

    def zpart(self, v) -> np.ndarray:
        return (self.gram @ np.array(v, dtype=np.int64)) % 4

    def embed_rows(self, pairs) -> np.ndarray:
        """pairs: list of (u_vec, v_vec) with u_vec, v_vec lists of GR elements (length n).
        Returns symplectic rows (a|b) over Z_4 of length 2N, N = m n."""
        rows = []
        for uvec, vvec in pairs:
            a = np.concatenate([self.xpart(u) for u in uvec])
            b = np.concatenate([self.zpart(v) for v in vvec])
            rows.append(np.concatenate([a, b]) % 4)
        return np.array(rows, dtype=np.int64)

    def trace_symplectic(self, uvec, vvec, uvec2, vvec2) -> int:
        s = self.zero
        for u, v, u2, v2 in zip(uvec, vvec, uvec2, vvec2):
            s = self.add(s, self.add(self.mul(u, v2), self.neg(self.mul(u2, v))))
        return self.trace(s)

    # ---- code builders
    def evaluation_code(self, points, band, twist):
        """Rows (ev(f) | ev(M f)) for f = x^j, j in band, M the twist (callable GR-vector -> GR-vector);
        M self-adjoint w.r.t. the trace form => isotropic.  Returns (rows, N)."""
        pairs = []
        for j in band:
            evf = [self.pow(p, j) for p in points]
            pairs.append((evf, twist(evf)))
        return self.embed_rows(pairs), self.m * len(points)

    def mult_twist(self, mu):
        return lambda evf: [self.mul(mu, u) for u in evf]

    def graph_code(self, C0_basis, Mmat):
        """Graph code {(c | M c) : c in C0}: C0_basis rows of GR elements (length n), Mmat a symmetric
        n x n GR matrix (self-adjoint for the trace form).  Returns (rows, N)."""
        n = len(C0_basis[0])
        pairs = []
        for c in C0_basis:
            Mc = []
            for i in range(n):
                acc = self.zero
                for j in range(n):
                    acc = self.add(acc, self.mul(Mmat[i][j], c[j]))
                Mc.append(acc)
            pairs.append((list(c), Mc))
        return self.embed_rows(pairs), self.m * n


@lru_cache(maxsize=None)
def gr(m: int) -> GR:
    return GR(m)


def random_symmetric_matrix(R: GR, n: int, rng) -> list[list[tuple]]:
    M = [[R.zero] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            e = R.elems[rng.randrange(len(R.elems))]
            M[i][j] = e; M[j][i] = e
    return M


def random_C0(R: GR, n: int, r: int, rng) -> list[list[tuple]]:
    return [[R.elems[rng.randrange(len(R.elems))] for _ in range(n)] for _ in range(r)]
