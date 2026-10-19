"""Exact solver for ALL transversal diagonal layers a code hosts (non-uniform per-qudit tables), and for
all across-block 2-body table layers.

Unknowns: per-qudit tables q_i : Z_d -> Z_M (M = 4d) with q_i(0) = 0, i.e. (d-1) n unknowns.  For a
layer Phi(x) = sum_i q_i(x_i):
  * hosting  <=>  Delta_a Phi(x) = sum_i [q_i(x_i + a_i) - q_i(x_i)] == 0 (mod M) for all X-generators a
                  and all points x of the coset spanning sets (support mode, exact for CSS; coset mode for
                  non-CSS with a free, soft-free logical quotient) — LINEAR in the unknowns;
  * level <= 3  <=>  third differences of each q_i constant — LINEAR.
The solution set is a Z_M-module computed with the Howell kernel.  The induced logical F(p) = Phi(x_p)
is linear in the unknowns as well, so escape of the logical gate (a non-vanishing condition) is decided
on a generating set / by sampling the module.  Likewise for pair tables T_(i,j) : Z_d^2 -> Z_M.

This generalizes the Phase-5 cubic solver (mod-4 polynomial Phi) to the full level-3 single-qudit table
space mod 16 and answers, for ANY code, "is there SOME transversal diagonal gate it hosts, and does any
hosted gate transmit universality?".
"""

from __future__ import annotations

import itertools as it
import math
import time

import numpy as np

from . import phase_tables as PT
from . import rgt
from .a2 import escapes_all_sectors
from .phase_tables import modulus
from .ringlinalg import howell_form, right_kernel_howell, module_size
from .symplectic import Z4Code, blocks


# --------------------------------------------------------------------------- unknown layout
class SingleQuditUnknowns:
    """q_i(v) for i in range(n), v in 1..d-1  ->  column index."""

    def __init__(self, n: int, d: int):
        self.n, self.d = n, d
        self.ncols = n * (d - 1)

    def col(self, i: int, v: int) -> int | None:
        v %= self.d
        return None if v == 0 else i * (self.d - 1) + (v - 1)

    def phi_row(self, x) -> np.ndarray:
        """Linear functional of Phi(x) = sum_i q_i(x_i)."""
        row = np.zeros(self.ncols, dtype=np.int64)
        for i, xi in enumerate(x):
            c = self.col(i, int(xi))
            if c is not None:
                row[c] += 1
        return row

    def tables(self, sol: np.ndarray) -> list[np.ndarray]:
        M = modulus(self.d)
        out = []
        for i in range(self.n):
            q = np.zeros(self.d, dtype=np.int64)
            q[1:] = sol[i * (self.d - 1):(i + 1) * (self.d - 1)] % M
            out.append(q)
        return out


def _coset_points(code: Z4Code, mode: str):
    """(base, gens) pairs whose spanning points carry the hosting conditions."""
    n = code.n
    if mode == "support":
        return [(np.zeros(n, dtype=np.int64), code.support_gens)]
    return [(x, code.x_gens) for x in rgt.logical_rep_points(code)]


def hosting_system(code: Z4Code, mode: str = "support", level3: bool = True) -> tuple[np.ndarray, SingleQuditUnknowns]:
    """Rows A such that A q == 0 (mod M) <=> the per-qudit layer q hosts (and is level <= 3)."""
    n, d = code.n, code.d
    M = modulus(d)
    U = SingleQuditUnknowns(n, d)
    rows = []
    for a in code.x_gens:
        a = np.asarray(a, dtype=np.int64)
        for base, gens in _coset_points(code, mode):
            pts = rgt.spanning_points(base, gens, d) if len(gens) else base[None, :] % d
            for x in pts:
                rows.append((U.phi_row((x + a) % d) - U.phi_row(x)) % M)
    if level3:
        for i in range(n):
            # third differences constant: D3 q_i(x) - D3 q_i(0) = 0 for x = 1..d-1
            def d3(x):
                r = np.zeros(U.ncols, dtype=np.int64)
                for mask in range(8):
                    bits = bin(mask).count("1")
                    y = (x + bits) % d
                    c = U.col(i, y)
                    if c is not None:
                        r[c] += (-1) ** (3 - bits)
                return r
            base = d3(0)
            for x in range(1, d):
                rows.append((d3(x) - base) % M)
    A = np.array(rows, dtype=np.int64).reshape(-1, U.ncols) if rows else np.zeros((0, U.ncols), dtype=np.int64)
    return A, U


def hosted_layer_module(code: Z4Code, mode: str = "support", level3: bool = True):
    """Generators (Howell) of the module of hosted per-qudit layers; returns (gens, U)."""
    A, U = hosting_system(code, mode, level3)
    M = modulus(code.d)
    if A.shape[0] == 0:
        gens = np.eye(U.ncols, dtype=np.int64)
    else:
        K = right_kernel_howell(A % M, M)
        gens = np.array(K, dtype=np.int64).reshape(-1, U.ncols) if len(K) else np.zeros((0, U.ncols), dtype=np.int64)
    return gens, U


def strict_level3_filter(gens: np.ndarray, U: SingleQuditUnknowns) -> list[np.ndarray]:
    """Module elements (generators + random combos) that are strict level-3 somewhere."""
    return [g for g in gens if any(PT.level(q, U.d) == 3 for q in U.tables(g))]


def analyze_code(code: Z4Code, mode: str | None = None, n_samples: int = 400, seed: int = 0, verbose: bool = False) -> dict:
    """Does the code host ANY transversal diagonal (level<=3) layer, and does any hosted layer induce an
    ESCAPING logical gate?  Exact module computation + sampling for the escape (non-linear) condition."""
    t0 = time.time()
    d = code.d
    M = modulus(d)
    if mode is None:
        mode = "support" if code.is_css else "cosets"
    exact = (code.is_css or (code.free_quotient_ok and code.soft == 0))
    gens, U = hosted_layer_module(code, mode)
    size = module_size(gens, M) if len(gens) else 1
    # Clifford-only layers: all q_i Clifford (level<=2) -> compute the size of that submodule for reference
    rep = {"label": code.label, "n": code.n, "k": code.k, "k_genuine": code.k_genuine, "soft": code.soft,
           "mode": mode, "exact": exact, "hosted_module_rank": int(len(gens)), "hosted_module_size": int(size)}
    rng = np.random.default_rng(seed)
    strict = 0; esc_phys = 0; esc_logical = 0; example = None
    cand = list(gens)
    for _ in range(n_samples):
        if len(gens) == 0:
            break
        c = rng.integers(0, M, size=len(gens))
        cand.append((c @ gens) % M)
    for sol in cand:
        tabs = U.tables(sol)
        if not any(PT.level(q, d) == 3 for q in tabs):
            continue
        strict += 1
        layer = rgt.TransversalLayer.per_qudit(tabs, d, "solved")
        if layer.escapes_physically():
            esc_phys += 1
            if 1 <= code.k_genuine <= 6:
                F = rgt.induced_logical(code, layer)
                if escapes_all_sectors(F, d):
                    esc_logical += 1
                    if example is None:
                        example = {"tables": [q.tolist() for q in tabs], "logical": F.reshape(-1).tolist(),
                                   "logical_level": PT.level(F, d)}
    rep.update({"candidates_tested": len(cand), "strict_level3_layers": strict, "physically_escaping": esc_phys,
                "escaping_logical": esc_logical, "example": example, "seconds": round(time.time() - t0, 1)})
    if verbose:
        print({k: v for k, v in rep.items() if k != "example"}, flush=True)
    return rep


# --------------------------------------------------------------------------- 2-body across-block tables
class PairUnknowns:
    """T(u,v) for u,v in Z_d with T(0,.) = T(.,0) = 0  ->  (d-1)^2 unknowns per pair; all pairs share
    one table (uniform across-block) or have independent tables."""

    def __init__(self, npairs: int, d: int, shared: bool = True):
        self.npairs, self.d, self.shared = npairs, d, shared
        self.per = (d - 1) ** 2
        self.ncols = self.per if shared else npairs * self.per

    def col(self, p: int, u: int, v: int):
        u %= self.d; v %= self.d
        if u == 0 or v == 0:
            return None
        base = 0 if self.shared else p * self.per
        return base + (u - 1) * (self.d - 1) + (v - 1)

    def phi_row(self, x, pairs):
        row = np.zeros(self.ncols, dtype=np.int64)
        for p, (i, j) in enumerate(pairs):
            c = self.col(p, int(x[i]), int(x[j]))
            if c is not None:
                row[c] += 1
        return row

    def table(self, sol, p=0):
        M = modulus(self.d)
        T = np.zeros((self.d, self.d), dtype=np.int64)
        base = 0 if self.shared else p * self.per
        for u in range(1, self.d):
            for v in range(1, self.d):
                T[u, v] = sol[base + (u - 1) * (self.d - 1) + (v - 1)] % M
        return T


def analyze_across_block_pairs(code: Z4Code, n_samples: int = 300, seed: int = 0, verbose: bool = False) -> dict:
    """All shared 2-body tables T on two copies of the code (qudit i of block 0 with qudit i of block 1)
    that host; report whether any is strict level-3 and induces an escaping logical."""
    t0 = time.time()
    d = code.d; M = modulus(d); n0 = code.n
    two = blocks(code, 2)
    pairs = [(i, n0 + i) for i in range(n0)]
    PU = PairUnknowns(n0, d, shared=True)
    mode = "support" if two.is_css else "cosets"
    rows = []
    for a in two.x_gens:
        a = np.asarray(a, dtype=np.int64)
        for base, gens in _coset_points(two, mode):
            pts = rgt.spanning_points(base, gens, d) if len(gens) else base[None, :] % d
            for x in pts:
                rows.append((PU.phi_row((x + a) % d, pairs) - PU.phi_row(x, pairs)) % M)
    # level <= 3 of the shared table: third unit differences constant (on the 2-qudit table)
    grid = list(it.product(range(d), repeat=2))
    def tab_row(u, v):
        r = np.zeros(PU.ncols, dtype=np.int64)
        c = PU.col(0, u, v)
        if c is not None:
            r[c] += 1
        return r
    for axes in it.combinations_with_replacement(range(2), 3):
        def d3(u, v):
            r = np.zeros(PU.ncols, dtype=np.int64)
            for mask in range(8):
                bits = bin(mask).count("1"); uu, vv = u, v
                for t, ax in enumerate(axes):
                    if mask >> t & 1:
                        if ax == 0: uu += 1
                        else: vv += 1
                r += (-1) ** (3 - bits) * tab_row(uu % d, vv % d)
            return r
        base = d3(0, 0)
        for (u, v) in grid[1:]:
            rows.append((d3(u, v) - base) % M)
    A = np.array(rows, dtype=np.int64).reshape(-1, PU.ncols)
    K = right_kernel_howell(A % M, M)
    gens = np.array(K, dtype=np.int64).reshape(-1, PU.ncols) if len(K) else np.zeros((0, PU.ncols), dtype=np.int64)
    rep = {"label": code.label, "n": code.n, "k": code.k, "pair_module_rank": int(len(gens)),
           "pair_module_size": int(module_size(gens, M)) if len(gens) else 1}
    rng = np.random.default_rng(seed)
    cand = list(gens) + [(rng.integers(0, M, size=len(gens)) @ gens) % M for _ in range(n_samples if len(gens) else 0)]
    strict = esc = escL = 0; example = None
    for sol in cand:
        T = PU.table(sol)
        if PT.level(T, d) != 3:
            continue
        strict += 1
        if escapes_all_sectors(T, d):
            esc += 1
            if 1 <= two.k_genuine <= 6:
                layer = rgt.TransversalLayer.across_block_pairs(n0, T, d)
                F = rgt.induced_logical(two, layer)
                if escapes_all_sectors(F, d):
                    escL += 1
                    if example is None:
                        example = {"table": T.tolist(), "logical_level": PT.level(F, d)}
    rep.update({"strict_level3": strict, "escaping_tables": esc, "escaping_logical": escL, "example": example,
                "seconds": round(time.time() - t0, 1)})
    if verbose:
        print({k: v for k, v in rep.items() if k != "example"}, flush=True)
    return rep
