"""Distance bounds and exact weight/distance computations.

This module provides four things (NOTES sec 7, Thm 4, sec 10 items 8/11):

  * ``delta_p``            -- the closed-form analytic distance Delta_p(m,r,w)
                             of a Manhattan-punctured Reed-Muller code (Thm 4).
  * ``min_distance``       -- exact minimum nonzero Hamming weight of an arbitrary
                             F_p linear code given a generator matrix.
  * ``weight_enumerator_at`` -- exact number of weight-d codewords of such a code.
  * ``A_d_logical_Z``      -- exact number of weight-d *logical*-Z operators of a
                             triorthogonal CSS code (used for the distillation
                             error model, NOTES eqs 38/39).

Counting strategy (low-weight, exact)
-------------------------------------
A weight-d codeword of a code C with parity-check matrix H is a set of d columns
of H carrying a linear dependency with *all-nonzero* coefficients:

    c has support T, |T| = d, all c_j != 0  <=>  sum_{j in T} c_j H[:,j] = 0.

So we iterate over the C(n,d) column d-subsets of H and, for each, enumerate the
full-support vectors in the (small) right null space of H[:, T].  This is exact
and is efficient whenever d is small and H has few rows -- i.e. exactly when the
code C is high-dimensional, which is the regime of the CSS distance/logical codes
(their parity checks are the small punctured/shortened RM generators).

Logical-Z operators (NOTES sec 3, sec 10 item 11)
-------------------------------------------------
For the CSS code built by ``triorthogonal.build_triorthogonal_code`` with
X-stabilizers ``X_stab = G0`` and Z-stabilizers ``Z_stab = Gp^perp``, with the
nesting G0 subset Gp subset G0^perp, the (nontrivial) **Z-logical operators** are
the codewords that commute with every X-stabilizer but are not Z-stabilizers:

    Z-logicals = G0^perp \\ Gp^perp     (commute with G0;  modulo Gp^perp).

Hence A_d counts the weight-d codewords of G0^perp that do NOT lie in Gp^perp.
(G0^perp has parity check G0 = X_stab; membership in Gp^perp is the test
Gp @ c^T == 0.)  This reproduces the A_d column of Table 3 exactly.
"""
from __future__ import annotations

import math
from itertools import combinations, product

import numpy as np

from .field import GFp
from .pnomial import pnomial_gt
from .triorthogonal import dual_matrix

# Budget guards (kept generous; the required oracle codes sit well inside them).
_DIRECT_BUDGET = 200_000          # enumerate p**k codewords directly below this
_COMB_BUDGET = 200_000_000        # refuse C(n,d) subset scans above this


# ---------------------------------------------------------------------------
# Analytic distance: Theorem 4
# ---------------------------------------------------------------------------
def delta_p(m, r, w, p) -> int:
    """Closed-form distance Delta_p(m,r,w) of a Manhattan-punctured RM code (NOTES Thm 4).

    With the degree decomposition r = alpha*(p-1) + beta, alpha = floor(r/(p-1)),
    beta = r mod (p-1), Theorem 4 (eqs near line 257) gives the unified form

        Delta_p(m,r,w) = sum_{j=0}^{p-beta-1} [ m-alpha-1, > (w-j) ]_p .

    The two closed branches are asserted as internal cross-checks:

        beta == 0:  Delta_p = [ m-alpha, > w ]_p
        beta != 0:  Delta_p = [ m-alpha, > w ]_p
                              - sum_{j=p-beta}^{p-1} [ m-alpha-1, > (w-j) ]_p

    both following from the generalized Pascal rule (pnomial Lemma 3).  Returns an
    exact python int; the formula is the code distance only when Delta_p > 0
    (otherwise the punctured generator is rank-deficient, NOTES Thm 4 note).
    """
    assert isinstance(p, int) and p >= 2, "p must be an integer prime >= 2"
    assert isinstance(m, int) and m >= 1, "m must be a positive integer"
    assert 0 <= r <= m * (p - 1), "r must lie in [0, m(p-1)]"

    alpha = r // (p - 1)
    beta = r - (p - 1) * alpha
    assert alpha <= m - 1, "Thm 4 requires alpha <= m-1 (so m-alpha-1 >= 0)"

    result = sum(pnomial_gt(m - alpha - 1, w - j, p) for j in range(0, p - beta))

    # Closed-form cross-checks (NOTES Thm 4, two branches).
    if beta == 0:
        assert result == pnomial_gt(m - alpha, w, p), "beta==0 closed form mismatch"
    else:
        correction = sum(
            pnomial_gt(m - alpha - 1, w - j, p) for j in range(p - beta, p)
        )
        assert result == pnomial_gt(m - alpha, w, p) - correction, (
            "beta!=0 closed form mismatch"
        )
    return int(result)


# ---------------------------------------------------------------------------
# Internal F_p linear algebra (exact, python big-int safe)
# ---------------------------------------------------------------------------
def _as_int_matrix(generator, p) -> np.ndarray:
    """Coerce a generator (galois or array-like) to a 2D int64 ndarray mod p."""
    M = np.asarray(generator)
    if M.ndim != 2:
        raise ValueError("generator must be a 2D matrix")
    return (M.astype(np.int64)) % p


def _right_null_space(cols, p):
    """Right null space of a small matrix given column-major as ``cols`` (r x d).

    ``cols`` is a list of ``d`` columns, each a length-r list of ints.  Returns a
    list of basis vectors (each length d) for { lambda in F_p^d : sum_j lambda_j
    cols[j] = 0 }.  Plain Gaussian elimination over F_p (p prime => inverses exist).
    """
    d = len(cols)
    r = len(cols[0]) if d else 0
    # Work on the r x d matrix A (row i, col j) = cols[j][i].
    A = [[cols[j][i] % p for j in range(d)] for i in range(r)]
    pivot_col_of_row = []
    pivot_cols = {}
    row = 0
    for col in range(d):
        piv = None
        for rr in range(row, r):
            if A[rr][col] % p != 0:
                piv = rr
                break
        if piv is None:
            continue
        A[row], A[piv] = A[piv], A[row]
        inv = pow(A[row][col], p - 2, p)
        A[row] = [(x * inv) % p for x in A[row]]
        for rr in range(r):
            if rr != row and A[rr][col] % p != 0:
                f = A[rr][col]
                A[rr] = [(A[rr][k] - f * A[row][k]) % p for k in range(d)]
        pivot_cols[col] = row
        pivot_col_of_row.append(col)
        row += 1
        if row == r:
            break
    free = [col for col in range(d) if col not in pivot_cols]
    basis = []
    for fcol in free:
        v = [0] * d
        v[fcol] = 1
        for pcol, prow in pivot_cols.items():
            v[pcol] = (-A[prow][fcol]) % p
        basis.append(v)
    return basis


def _full_support_kernel_vectors(basis, p, d):
    """Yield every all-nonzero length-d vector in the F_p span of ``basis``."""
    t = len(basis)
    if t == 0:
        return
    for coeffs in product(range(p), repeat=t):
        if not any(coeffs):
            continue
        vec = [0] * d
        for bi, cf in enumerate(coeffs):
            if cf:
                row = basis[bi]
                for k in range(d):
                    vec[k] = (vec[k] + cf * row[k]) % p
        if all(x != 0 for x in vec):
            yield vec


def _iter_weight_d_from_parity(H, p, d, collect, early=False):
    """Count / collect weight-d codewords of the code with parity check ``H``.

    ``H`` is an int ndarray (r x n).  A weight-d codeword is a full-support
    dependency among d columns of H (see module docstring).  If ``collect`` is
    True, return a list of length-n codeword vectors (ndarray); otherwise return
    the integer count.  ``early`` (count mode only) returns 1 as soon as one is
    found (used by ``min_distance``).
    """
    r, n = H.shape
    if d < 0 or d > n:
        return [] if collect else 0
    if d == 0:
        return [np.zeros(n, dtype=np.int64)] if collect else 1

    # Degenerate code = full space F_p^n (no parity checks): every full-support
    # d-subset is a codeword.  Count = C(n,d) * (p-1)^d.
    if r == 0:
        if early and not collect:
            return 1
        if not collect:
            return math.comb(n, d) * (p - 1) ** d
        out = []
        for T in combinations(range(n), d):
            for vals in product(range(1, p), repeat=d):
                c = np.zeros(n, dtype=np.int64)
                for k, col in enumerate(T):
                    c[col] = vals[k]
                out.append(c)
        return out

    Hcols = [H[:, j].tolist() for j in range(n)]
    count = 0
    out = []
    for T in combinations(range(n), d):
        sub = [Hcols[j] for j in T]
        basis = _right_null_space(sub, p)
        if not basis:
            continue
        for lam in _full_support_kernel_vectors(basis, p, d):
            if collect:
                c = np.zeros(n, dtype=np.int64)
                for k, col in enumerate(T):
                    c[col] = lam[k]
                out.append(c)
            else:
                count += 1
                if early:
                    return 1
    return out if collect else count


def _code_basis(generator, p):
    """Return (basis ndarray (k x n) of independent rows, k, n) for the row span."""
    GF = GFp(p)
    G = GF(_as_int_matrix(generator, p))
    R = np.asarray(G.row_reduce()).astype(np.int64) % p
    nz = [row for row in R if np.any(row != 0)]
    n = G.shape[1]
    if not nz:
        return np.zeros((0, n), dtype=np.int64), 0, n
    B = np.array(nz, dtype=np.int64)
    return B, B.shape[0], n


# ---------------------------------------------------------------------------
# Exact minimum distance
# ---------------------------------------------------------------------------
def min_distance(generator, p, max_weight=None) -> int:
    """Exact minimum nonzero Hamming weight of the F_p code spanned by ``generator``.

    Two exact strategies (NOTES sec 10 item 8):

      * direct enumeration of all p**k codewords when k = dim is small
        (p**k <= internal budget);
      * otherwise an increasing-weight column-dependency scan on the parity
        check (efficient when the code is high-dimensional, e.g. a punctured-RM
        dual), returning the first weight d that admits a codeword.

    ``max_weight`` bounds the dependency scan; if no codeword of weight
    <= max_weight exists the true distance is *not* under-reported -- a
    ValueError is raised instead.  Raises ValueError for the trivial (zero) code.
    """
    assert isinstance(p, int) and p >= 2, "p must be an integer prime >= 2"
    B, k, n = _code_basis(generator, p)
    if k == 0:
        raise ValueError("trivial (zero-dimensional) code has no nonzero codeword")

    if p ** k <= _DIRECT_BUDGET:
        best = None
        for coeffs in product(range(p), repeat=k):
            if not any(coeffs):
                continue
            word = np.zeros(n, dtype=np.int64)
            for i, cf in enumerate(coeffs):
                if cf:
                    word = (word + cf * B[i]) % p
            w = int(np.count_nonzero(word))
            if best is None or w < best:
                best = w
                if best == 1:
                    break
        return best

    # High-dimensional code: scan increasing weight via the parity check.
    H = _as_int_matrix(dual_matrix(generator, p), p)
    limit = n if max_weight is None else min(max_weight, n)
    for d in range(1, limit + 1):
        if _iter_weight_d_from_parity(H, p, d, collect=False, early=True):
            return d
    if max_weight is not None:
        raise ValueError(
            f"minimum distance exceeds max_weight={max_weight}; not certified"
        )
    raise ValueError("no nonzero codeword found (unexpected for a nontrivial code)")


# ---------------------------------------------------------------------------
# Weight enumerator at a fixed weight
# ---------------------------------------------------------------------------
def weight_enumerator_at(generator, p, d) -> int:
    """Exact number of weight-d codewords of the F_p code spanned by ``generator``.

    Uses direct codeword enumeration when the dimension is small, otherwise the
    parity-check column-dependency count (NOTES sec 7).  ``d == 0`` returns 1
    (the zero codeword).  Raises NotImplementedError if the C(n,d) subset scan
    would exceed the internal budget.
    """
    assert isinstance(p, int) and p >= 2, "p must be an integer prime >= 2"
    B, k, n = _code_basis(generator, p)
    if d == 0:
        return 1
    if d < 0 or d > n:
        return 0

    if p ** k <= _DIRECT_BUDGET:
        count = 0
        for coeffs in product(range(p), repeat=k):
            if not any(coeffs):
                continue
            word = np.zeros(n, dtype=np.int64)
            for i, cf in enumerate(coeffs):
                if cf:
                    word = (word + cf * B[i]) % p
            if int(np.count_nonzero(word)) == d:
                count += 1
        return count

    if math.comb(n, d) > _COMB_BUDGET:
        raise NotImplementedError(
            f"weight enumerator: C({n},{d}) = {math.comb(n, d)} exceeds budget "
            f"{_COMB_BUDGET}; exact counting infeasible for this code"
        )
    H = _as_int_matrix(dual_matrix(generator, p), p)
    return int(_iter_weight_d_from_parity(H, p, d, collect=False))


# ---------------------------------------------------------------------------
# Weight-d logical-Z operator count
# ---------------------------------------------------------------------------
def A_d_logical_Z(code, p, d) -> int:
    """Number of weight-d logical-Z operators of a triorthogonal CSS code (NOTES sec 3, 11).

    ``code`` is the dict returned by ``triorthogonal.build_triorthogonal_code``
    (keys ``X_stab`` = G0, ``Gp``).  The Z-logical operators are the codewords of
    G0^perp (commute with all X-stabilizers G0) that are NOT Z-stabilizers, i.e.
    not in Gp^perp (NOTES sec 3, nesting G0 subset Gp subset G0^perp):

        A_d = #{ weight-d c in G0^perp : Gp @ c^T != 0 }.

    G0^perp has parity check G0, so weight-d codewords of G0^perp are enumerated
    by the column-dependency method on G0; each is then tested against Gp.  This
    reproduces the A_d column of Table 3 (NOTES line ~381).

    Raises NotImplementedError when the exact C(n,d) subset scan would exceed the
    internal budget (the test layer skips-with-reason in that case); it never
    returns a guessed value.
    """
    assert isinstance(p, int) and p >= 2, "p must be an integer prime >= 2"
    X = _as_int_matrix(code["X_stab"], p)   # G0  (parity check of G0^perp)
    Gp = _as_int_matrix(code["Gp"], p)      # punctured generator
    n = X.shape[1]
    assert Gp.shape[1] == n, "X_stab and Gp must share the block length n"

    if d < 0 or d > n:
        return 0
    if math.comb(n, d) > _COMB_BUDGET:
        raise NotImplementedError(
            f"A_d: C({n},{d}) = {math.comb(n, d)} exceeds budget {_COMB_BUDGET}; "
            f"exact weight-d logical enumeration infeasible for this code"
        )

    words = _iter_weight_d_from_parity(X, p, d, collect=True)
    if not words:
        return 0
    W = np.array(words, dtype=np.int64)              # (#codewords, n)
    syndromes = (W @ Gp.T) % p                        # (#codewords, rows of Gp)
    not_stabilizer = np.any(syndromes != 0, axis=1)   # Gp c^T != 0  => true logical
    return int(np.count_nonzero(not_stabilizer))
