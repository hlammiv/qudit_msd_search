"""Weakly-self-dual CSS codes over the cyclic ring Z_d (d = 2^k), with the code-level anti-collapse
test.

A CSS code is fixed by its X-stabilizer module M_X (rows of a generator matrix over Z_d).  We take
the weakly-self-dual choice M_X subseteq M_X^perp with Z-stabilizers = M_X, so the logical degrees of
freedom live in the quotient  M_X^perp / M_X,  an abelian group of order d^k.

The CODE-LEVEL anti-collapse certificate (the part that distinguishes a genuine single-ring-qudit
code from "two qubits in disguise"):  for k = 1 the logical group must be the **cyclic** group Z_d,
not Z_2 x Z_{d/2}.  If M_X^perp / M_X factorizes into smaller cyclic groups, the "logical qudit" is a
product of smaller qudits and the construction has collapsed.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .ringlinalg import enumerate_module, dual_module, module_size


@dataclass
class CSSCode:
    n: int
    d: int
    gens: tuple                      # X-stabilizer generators (rows over Z_d)
    mx: list = field(repr=False)     # elements of M_X (X-stabilizer module)
    mxperp: list = field(repr=False)  # elements of M_Z^perp = valid logical-rep space (== M_X^perp when M_Z=M_X)
    k: int = 0                       # number of logical qudits
    cyclic: bool = False             # logical group is cyclic Z_d (k==1) -> code anti-collapse OK
    logical_x: tuple | None = None   # coset rep generating the cyclic logical group
    cosets: list = field(default_factory=list, repr=False)  # canonical reps r_0..r_{d-1} (k==1, cyclic)
    mz: list | None = field(default=None, repr=False)  # M_Z (Z-stabilizer module); None for weakly-self-dual
    distance: int = 0                # min(d_X, d_Z); set by build_css_general


def is_self_orthogonal(gens, d: int) -> bool:
    """M_X subseteq M_X^perp, i.e. a . a' == 0 (mod d) for all generator pairs (incl. a == a')."""
    G = np.atleast_2d(np.array(gens, dtype=object)) % d
    return bool(np.all((G @ G.T) % d == 0))


def build_css(gens, d: int) -> CSSCode | None:
    """Build the weakly-self-dual CSS code from X-stabilizer generators, or None if not self-orthogonal."""
    G = np.atleast_2d(np.array(gens, dtype=object)) % d
    n = G.shape[1]
    if not is_self_orthogonal(G, d):
        return None
    mx = enumerate_module(G, d)
    mx_set = set(mx)
    mxperp = dual_module(G, d)
    sz_mx = len(mx)
    sz_perp = len(mxperp)
    # number of logical qudits: d^k = [M_X^perp : M_X]
    ratio = sz_perp // sz_mx
    k = 0
    r = ratio
    while r > 1 and r % d == 0:
        r //= d
        k += 1
    if r != 1:
        # ratio is not a pure power of d -> mixed structure; record k via log but flag non-standard
        k = -1

    code = CSSCode(n=n, d=d, gens=tuple(map(tuple, G.tolist())), mx=mx, mxperp=mxperp, k=k)

    if k == 1:
        # find a coset rep g of order d in the quotient M_X^perp / M_X  (=> cyclic Z_d)
        for g in mxperp:
            if g in mx_set:
                continue
            order = None
            for t in range(1, d + 1):
                tg = tuple((t * np.array(g)) % d)
                if tg in mx_set:
                    order = t
                    break
            if order == d:
                code.cyclic = True
                code.logical_x = tuple(g)
                code.cosets = [tuple((ell * np.array(g)) % d) for ell in range(d)]
                break
    return code


def _cyclic_cosets(mx_set, mzperp, d):
    """If mzperp/mx is cyclic Z_d (order-d element exists), return (logical_x, [r_0..r_{d-1}]); else None."""
    for g in mzperp:
        if g in mx_set:
            continue
        order = next((t for t in range(1, d + 1) if tuple((t * np.array(g)) % d) in mx_set), None)
        if order == d:
            return tuple(g), [tuple((ell * np.array(g)) % d) for ell in range(d)]
    return None


def _min_weight_outside(elems, sub_set, d):
    """Min Hamming weight over elems not in sub_set (a logical-operator distance)."""
    best = 10 ** 9
    for y in elems:
        if y in sub_set:
            continue
        w = sum(1 for c in y if c % d != 0)
        if 0 < w < best:
            best = w
    return best


def build_css_general(A, B, d: int, mxperp_pre=None, k1_only: bool = True) -> CSSCode | None:
    """General CSS code over Z_d from X-stabilizers A and Z-stabilizers B (M_X subseteq M_Z^perp).

    Logical reps live in M_Z^perp / M_X.  Distance = min(d_X, d_Z) with d_X = min wt in M_Z^perp \\ M_X,
    d_Z = min wt in M_X^perp \\ M_Z.  Returns None unless the commutation A B^T == 0 holds.

    Fast path (``k1_only``): since |M||M^perp| = d^n for the (perfect) standard pairing, k=1 is exactly
    |M_X|*|M_Z| == d^{n-1}.  We check this from the cheap module sizes and skip the expensive d^n dual
    enumeration entirely when the code cannot be k=1.  ``mxperp_pre`` lets the caller pass M_X^perp
    (already computed to sample B from it) so it is not recomputed.
    """
    A = np.atleast_2d(np.array(A, dtype=object)) % d
    B = np.atleast_2d(np.array(B, dtype=object)) % d
    n = A.shape[1]
    if not np.all((A @ B.T) % d == 0):       # X/Z commutation
        return None
    mx = enumerate_module(A, d)
    mz = enumerate_module(B, d)
    if k1_only and len(mx) * len(mz) != d ** (n - 1):   # cheap k=1 prune (before any d^n dual)
        return None
    mx_set, mz_set = set(mx), set(mz)
    mzperp = dual_module(B, d)               # valid logical-rep space
    ratio = len(mzperp) // len(mx)
    k, r = 0, ratio
    while r > 1 and r % d == 0:
        r //= d
        k += 1
    if r != 1:
        k = -1
    code = CSSCode(n=n, d=d, gens=tuple(map(tuple, A.tolist())), mx=mx, mxperp=mzperp, k=k, mz=mz)
    if k == 1:
        cc = _cyclic_cosets(mx_set, mzperp, d)
        if cc is not None:
            code.cyclic = True
            code.logical_x, code.cosets = cc
            mxperp = mxperp_pre if mxperp_pre is not None else dual_module(A, d)
            dX = _min_weight_outside(mzperp, mx_set, d)
            dZ = _min_weight_outside(mxperp, mz_set, d)
            code.distance = min(dX, dZ)
    return code


def code_distance(code: CSSCode) -> int:
    """Distance = min Hamming weight of a nontrivial logical operator (weakly-self-dual: M_X^perp \\ M_X).

    Distance 1 means a logical operator acts on a single physical qudit (the logical qudit is
    effectively unencoded) — such a code cannot distill.  A genuine distillation code needs d >= 2.
    """
    mx_set = set(code.mx)
    best = code.n + 1
    for y in code.mxperp:
        if y in mx_set:
            continue
        w = sum(1 for c in y if c % code.d != 0)
        if 0 < w < best:
            best = w
    return best


def logical_group_invariants(code: CSSCode) -> list[int]:
    """Element orders present in M_X^perp / M_X (diagnostic for the cyclic-vs-product question)."""
    mx_set = set(code.mx)
    orders = set()
    for g in code.mxperp:
        for t in range(1, code.d + 1):
            if tuple((t * np.array(g)) % code.d) in mx_set:
                orders.add(t)
                break
    return sorted(orders)
