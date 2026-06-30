"""Structured A_d enumerator for punctured Reed-Muller / triorthogonal qudit codes.

For a punctured triorthogonal code the quantum code is ``G0^perp`` = the *punctured*
generalized Reed-Muller code ``RM_p(rtilde, m)`` restricted to the surviving
coordinates, where ``S`` (``|S| = k``) is the set of punctured points
(``rtilde = m(p-1) - r - 1``).  Whenever the punctured-column submatrix is full rank
the restriction map ``RM_p(rtilde,m) -> G0^perp`` is a *bijection* (its kernel
``{c : supp(c) subset S}`` is trivial), so for any nonzero ``c in RM_p(rtilde,m)`` the
punctured weight is ``|supp(c) \\ S|`` and

    d   = min_{c != 0} |supp(c) \\ S|
    A_d = #{ c in RM_p(rtilde,m), c != 0 : |supp(c) \\ S| = d }.

Naively this needs the MacWilliams engine (cost ``p^{dim G0} = p^{D-k}``) or a direct
``C(n,d)`` column scan -- both blow up for the high-distance / low-redundancy codes.

The structured method (this module)
-----------------------------------
Every weight-d codeword of ``G0^perp`` lifts to a unique low-weight ``RM_p(rtilde,m)``
codeword ``c``, and such codewords are *geometrically structured*: ``c`` is supported on
an affine flat ``F = aff.span(supp c)``, and restricted to ``F ~= F_p^j`` it is a codeword
of the (much smaller) code

    C_F = RM_p( rtilde - (m-j)(p-1),  j )            [the codewords of RM_p(rtilde,m)
                                                      supported inside the j-flat F].

We therefore enumerate weight-d codewords by their MINIMAL affine span dimension ``j``:
for each ``j``-flat ``F`` we find the codewords of ``C_F`` with exactly ``d`` non-zeros on
the surviving points ``F \\ S`` (a meet-in-the-middle column-dependency search on the small
parity check of ``C_F`` restricted to ``F``), and keep those whose support spans ``F``
exactly (``aff.dim == j``) -- the minimal-flat condition makes every codeword counted once.

Cost scales with the number of low-dimensional flats and the size of the (tiny) restricted
codes, NOT with ``p^{D-k}`` or ``C(n,d)``.  Summed over all ``j = 2 .. m`` the count is
EXACT; capping ``jmax < m`` gives the contribution from codewords of affine span ``<= jmax``
(a certified lower bound on ``A_d``, exact when no higher-span codeword contributes).

This is the algebraic structure behind the qutrit ``A_d`` values of arXiv:2510.10852.
"""
from __future__ import annotations

import itertools
from collections import Counter
from functools import lru_cache

import numpy as np

from .reedmuller import points, rm_generator
from .triorthogonal import build_triorthogonal_code

__all__ = ["structured_ad"]


# ---------------------------------------------------------------------------
# pure-numpy exact F_p linear algebra (no galois per-call overhead)
# ---------------------------------------------------------------------------
def _rref(M, p):
    """Reduced row-echelon form over F_p (in place on a copy). Returns (R, pivot_cols)."""
    R = (np.asarray(M, dtype=np.int64) % p).copy()
    rows, cols = R.shape
    pivots = []
    pr = 0
    for c in range(cols):
        sel = -1
        for rr in range(pr, rows):
            if R[rr, c] % p != 0:
                sel = rr
                break
        if sel < 0:
            continue
        R[[pr, sel]] = R[[sel, pr]]
        inv = pow(int(R[pr, c]), p - 2, p)
        R[pr] = (R[pr] * inv) % p
        for rr in range(rows):
            if rr != pr and R[rr, c] % p != 0:
                R[rr] = (R[rr] - R[rr, c] * R[pr]) % p
        pivots.append(c)
        pr += 1
        if pr == rows:
            break
    return R, pivots


def _rank(M, p):
    """Exact rank of M over F_p."""
    _, pivots = _rref(M, p)
    return len(pivots)


def _null_space(G, p):
    """Right null space basis of G (k x n) over F_p, as rows of an (n-rank) x n array."""
    G = np.asarray(G, dtype=np.int64) % p
    k, n = G.shape
    R, pivots = _rref(G, p)
    pivset = set(pivots)
    free = [c for c in range(n) if c not in pivset]
    basis = np.zeros((len(free), n), dtype=np.int64)
    for i, f in enumerate(free):
        basis[i, f] = 1
        for prow, pc in enumerate(pivots):
            basis[i, pc] = (-R[prow, f]) % p
    return basis


def _inv_modp(B, p):
    """Inverse of a square full-rank matrix B over F_p via Gauss-Jordan on [B | I]."""
    B = np.asarray(B, dtype=np.int64) % p
    n = B.shape[0]
    M = np.concatenate([B, np.eye(n, dtype=np.int64)], axis=1) % p
    for c in range(n):
        sel = -1
        for rr in range(c, n):
            if M[rr, c] % p != 0:
                sel = rr
                break
        if sel < 0:
            raise ValueError("singular matrix in _inv_modp")
        M[[c, sel]] = M[[sel, c]]
        inv = pow(int(M[c, c]), p - 2, p)
        M[c] = (M[c] * inv) % p
        for rr in range(n):
            if rr != c and M[rr, c] % p != 0:
                M[rr] = (M[rr] - M[rr, c] * M[c]) % p
    return M[:, n:] % p


def _systematic_lift(Gsurv, Gj, p):
    """Per-flat systematic decoder: codeword-values-on-survivors -> full codeword on the flat.

    ``Gsurv`` (dimCF x |surv|) is the flat generator restricted to surviving columns.  Because
    the global puncture submatrix is full rank, the restriction ``C_F -> C_F|surv`` is injective,
    so ``Gsurv`` has full row rank dimCF and an information set exists *among the survivors*.
    Returns ``(pivots, T)``: ``pivots`` = the dimCF survivor-local column indices of an info set,
    and ``T = B^{-1} @ Gj`` (dimCF x p^j) with ``B = Gsurv[:, pivots]``.  For any codeword the full
    evaluation on the whole flat is ``c[pivots] @ T`` -- one rref+inverse per flat replaces the
    per-codeword linear solve (the old hot path)."""
    _, pivots = _rref(Gsurv, p)
    B = np.asarray(Gsurv, dtype=np.int64)[:, pivots] % p
    Binv = _inv_modp(B, p)
    T = (Binv @ (np.asarray(Gj, dtype=np.int64) % p)) % p
    return pivots, T


def _solve(A, b, p, ncols):
    """One solution u (length ncols) of A u = b over F_p (A: r x ncols), assumed consistent."""
    A = np.asarray(A, dtype=np.int64) % p
    b = np.asarray(b, dtype=np.int64) % p
    M = np.concatenate([A, b.reshape(-1, 1)], axis=1) % p
    rows = M.shape[0]
    pivots = []
    pr = 0
    for c in range(ncols):
        sel = -1
        for rr in range(pr, rows):
            if M[rr, c] % p != 0:
                sel = rr
                break
        if sel < 0:
            continue
        M[[pr, sel]] = M[[sel, pr]]
        inv = pow(int(M[pr, c]), p - 2, p)
        M[pr] = (M[pr] * inv) % p
        for rr in range(rows):
            if rr != pr and M[rr, c] % p != 0:
                M[rr] = (M[rr] - M[rr, c] * M[pr]) % p
        pivots.append(c)
        pr += 1
        if pr == rows:
            break
    u = np.zeros(ncols, dtype=np.int64)
    for i, c in enumerate(pivots):
        u[c] = M[i, -1] % p
    return u


# ---------------------------------------------------------------------------
# affine flats of F_p^m
# ---------------------------------------------------------------------------
def _subspaces(m, j, p):
    """Yield every j-dim linear subspace of F_p^m once as a (j x m) RREF basis."""
    for pivots in itertools.combinations(range(m), j):
        pivotset = set(pivots)
        free_per_row = [[c for c in range(pc + 1, m) if c not in pivotset] for pc in pivots]
        ranges = [itertools.product(range(p), repeat=len(fr)) for fr in free_per_row]
        for combo in itertools.product(*ranges):
            basis = np.zeros((j, m), dtype=np.int64)
            for i, pc in enumerate(pivots):
                basis[i, pc] = 1
                for c, v in zip(free_per_row[i], combo[i]):
                    basis[i, c] = v
            yield basis, pivots


def _flats(m, j, p):
    """Yield colmap for each j-flat: colmap[t] = global point-index of the F_p^j point t.

    ``t`` ranges over ``points(j, p)`` (canonical order, x_1 least significant), so colmap
    aligns the columns of ``rm_generator(rr, j, p)`` with the flat's global point indices.
    Cosets of each subspace are generated in one batched numpy step (fast for large m).
    """
    POWm = (p ** np.arange(m)).astype(np.int64)
    Fj = points(j, p)                                   # (p^j, j)
    for basis, pivots in _subspaces(m, j, p):
        sub = (Fj @ basis) % p                          # (p^j, m)
        nonpiv = [c for c in range(m) if c not in pivots]
        # all coset representatives: 0 on pivot columns, free on the rest
        nfree = len(nonpiv)
        reps = np.array(list(itertools.product(range(p), repeat=nfree)), dtype=np.int64)
        origins = np.zeros((reps.shape[0], m), dtype=np.int64)
        if nfree:
            origins[:, nonpiv] = reps
        flats = (sub[None, :, :] + origins[:, None, :]) % p   # (R, p^j, m)
        codes = (flats @ POWm).astype(np.int64)               # (R, p^j)
        for row in codes:
            yield row


# ---------------------------------------------------------------------------
# meet-in-the-middle: weight-d codewords of H^perp (H = parity check)
# ---------------------------------------------------------------------------
# Overflow-proof syndrome encoding.  The old int64 p-adic ``powers @ s`` is bijective but
# raises ``OverflowError`` once ``p**rows > 2**63`` -- which happens at m=4, j>=3 where the
# restricted code ``C_F`` has large redundancy ``rows`` (e.g. 5**29).  The MITM only needs
# syndrome EQUALITY (``left_sum + right_sum = 0``), so we encode each length-``rows`` syndrome
# vector with a 64-bit polynomial HASH (uint64, wraps mod 2**64) and resolve every hash hit by
# verifying the actual parity ``H @ codeword == 0`` (a hash collision is then simply rejected).
# This is strictly more general than the int64 path and gives identical results.
_HASH_BASE = np.uint64(0x9E3779B97F4A7C15)  # fixed 64-bit golden-ratio multiplier


@lru_cache(maxsize=None)
def _base_powers(rows):
    """``[B^0, B^1, ..., B^(rows-1)]`` as uint64 (mod 2**64); the polynomial-hash basis.

    Computed in Python big-ints reduced mod 2**64 so the modular reduction is explicit (no
    noisy numpy uint64-overflow RuntimeWarning); the resulting wrap is exactly mod 2**64.
    """
    B = int(_HASH_BASE)
    MASK = (1 << 64) - 1
    vals = [1]
    for _ in range(1, rows):
        vals.append((vals[-1] * B) & MASK)
    return np.array(vals, dtype=np.uint64)


def _hash_cols(s, basepows):
    """Hash each column of ``s`` (rows x NC, entries in [0,p)) -> (NC,) uint64 polynomial hash.

    The per-element products and the column sum are computed in uint64 and wrap mod 2**64
    (an intentional modular hash); errstate silences numpy's overflow warning for the wrap.
    """
    with np.errstate(over="ignore"):
        return (basepows[:, None] * s.astype(np.uint64)).sum(axis=0)


def _half_sums(H, basepows, d_half, p, negate=False):
    """For all d_half-subsets x all-nonzero-coeffs: hashed syndrome, support, coeffs.

    Returns (hashes (N,) uint64, supps (N,d_half), coefs (N,d_half)).  Vectorised with numpy.
    ``negate`` encodes the *negated* syndrome ``(-s) mod p`` (componentwise) BEFORE hashing --
    used for the right half so that a left/right hash match means (after verification)
    ``left_sum + right_sum = 0``.  Coefficients returned are the actual (positive) codeword coeffs.
    """
    rows, n = H.shape
    combos = np.fromiter(
        itertools.chain.from_iterable(itertools.combinations(range(n), d_half)),
        dtype=np.int64,
    )
    if combos.size == 0:
        return (np.empty(0, np.uint64), np.empty((0, d_half), np.int64),
                np.empty((0, d_half), np.int64))
    combos = combos.reshape(-1, d_half)                       # (NC, d_half)
    Hsel = H[:, combos]                                        # (rows, NC, d_half)
    code_chunks, supp_chunks, coef_chunks = [], [], []
    for cf in itertools.product(range(1, p), repeat=d_half):
        cfa = np.array(cf, dtype=np.int64)
        s = (Hsel * cfa[None, None, :]).sum(axis=2) % p        # (rows, NC)
        if negate:
            s = (p - s) % p                                    # encode -syndrome (componentwise)
        code_chunks.append(_hash_cols(s, basepows))            # (NC,) uint64 hash
        supp_chunks.append(combos)
        coef_chunks.append(np.broadcast_to(cfa, combos.shape))
    return (np.concatenate(code_chunks),
            np.concatenate(supp_chunks),
            np.concatenate(coef_chunks))


def _mitm(H, d, p):
    """All weight-d codewords of the code with parity check H (rows x n).

    Returns a list of (support_tuple, values_tuple).  Splits d = d1 + d2, hashes each
    F_p^rows partial syndrome to one uint64 (overflow-proof), matches left/right halves by
    sorted searchsorted, and VERIFIES every hash hit against the true parity ``H @ cw == 0``
    (rejecting the rare hash collision) -- fully vectorised in numpy.
    """
    rows, n = H.shape
    if d > n or n == 0 or rows == 0:
        return []
    H = np.asarray(H, dtype=np.int64) % p
    basepows = _base_powers(rows)

    d1, d2 = d // 2, d - d // 2
    Lcode, Lsupp, Lcoef = _half_sums(H, basepows, d1, p)
    if Lcode.size == 0:
        return []
    order = np.argsort(Lcode, kind="stable")
    Lcode, Lsupp, Lcoef = Lcode[order], Lsupp[order], Lcoef[order]

    target, Rsupp, Rcoef = _half_sums(H, basepows, d2, p, negate=True)  # hash(-right_sum)
    pos = np.searchsorted(Lcode, target, side="left")
    valid = pos < Lcode.size
    match = np.zeros(target.shape, dtype=bool)
    match[valid] = Lcode[pos[valid]] == target[valid]

    out = set()
    for ri in np.nonzero(match)[0]:
        t = target[ri]
        lo = np.searchsorted(Lcode, t, side="left")
        hi = np.searchsorted(Lcode, t, side="right")
        rsupp = Rsupp[ri]
        rset = set(int(x) for x in rsupp)
        for li in range(lo, hi):
            lsupp = Lsupp[li]
            if rset.isdisjoint(int(x) for x in lsupp):
                vals = {}
                for c, a in zip(lsupp, Lcoef[li]):
                    vals[int(c)] = int(a)
                for c, a in zip(rsupp, Rcoef[ri]):
                    vals[int(c)] = int(a)
                keys = tuple(sorted(vals))
                # verify the actual parity (reject a hash collision): H[:,keys] @ vals == 0
                cols = np.array(keys, dtype=np.int64)
                coef = np.array([vals[k] for k in keys], dtype=np.int64)
                if np.any((H[:, cols] @ coef) % p):
                    continue
                out.add((keys, tuple(vals[k] for k in keys)))
    return list(out)


def _full_space_weight_d(n, d, p, budget=5_000_000):
    """Weight-d full-support codewords when the restricted code is the whole space F_p^n.

    Happens when the survivors exactly span C_F (no parity checks); every weight-d
    pattern is then a codeword.  Returns (support_tuple, values_tuple) list.
    """
    from math import comb
    if comb(n, d) * (p - 1) ** d > budget:
        raise ValueError(
            f"full-space restricted code with C({n},{d})*(p-1)^{d} weight-d patterns "
            f"exceeds budget {budget}; cannot enumerate this flat exactly"
        )
    out = []
    for sup in itertools.combinations(range(n), d):
        for vals in itertools.product(range(1, p), repeat=d):
            out.append((sup, vals))
    return out


# ---------------------------------------------------------------------------
# minimum weight of a codeword whose support spans the whole flat (pruning bound)
# ---------------------------------------------------------------------------
@lru_cache(maxsize=None)
def _min_span_weight(rr, j, p):
    """Min Hamming weight of a codeword of RM_p(rr,j) with affine span = F_p^j (dim exactly j).

    A *sound* pruning bound: a span-`j` codeword contributing to A_d has
    `|supp ∩ S| = w − d ≥ min_span_weight(j) − d`, so a `j`-flat with fewer than that many
    puncture points cannot host one.  Computed EXACTLY by brute enumeration when the (tiny)
    restricted code is small; for the large ones it uses the closed form
    ``p·(j+1)`` for the ``p=3`` ``b=0`` family (``rr = 2j−4``: gives 9,12,15,18,… — validated
    against the Table-3 dim-histograms), falling back to the safe ``d_RM`` lower bound.
    """
    from .reedmuller import d_rm
    Gj = np.asarray(rm_generator(rr, j, p)).astype(np.int64) % p
    dim = Gj.shape[0]
    if p ** dim <= 200_000:
        pts = points(j, p)
        msgs = np.array(list(itertools.product(range(p), repeat=dim)), dtype=np.int64)
        cws = (msgs @ Gj) % p
        weights = np.count_nonzero(cws, axis=1)
        order = np.argsort(weights, kind="stable")
        for idx in order:
            w = int(weights[idx])
            if w == 0:
                continue
            supp = np.nonzero(cws[idx])[0]
            P = pts[supp]
            adim = _rank((P[1:] - P[0]) % p, p) if P.shape[0] > 1 else 0
            if adim == j:
                return w
        return j * (p - 1) + 1
    if p == 3 and rr == 2 * j - 4:
        return p * (j + 1)
    return d_rm(rr, j, p)


# ---------------------------------------------------------------------------
# distance
# ---------------------------------------------------------------------------
def _distance(G0, p):
    # G0 is the parity check of the code G0^perp whose distance we want, so the
    # distance is the minimum number of F_p-dependent columns of G0 (certified MITM).
    from .mindist import min_dependent_columns
    return int(min_dependent_columns(np.asarray(G0).astype(np.int64) % p, p))


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def structured_ad(p, m, r, puncture_columns, jmax=None, distance=None, return_codewords=False):
    """Structured A_d of the punctured RM / triorthogonal code (no MacWilliams, no C(n,d) scan).

    Parameters
    ----------
    p, m, r : ints
        Prime ``p``, number of RM variables ``m``, base RM degree ``r`` (= r_max for the
        paper's codes); the dual degree is ``rtilde = m(p-1) - r - 1``.
    puncture_columns : iterable of 1-indexed column indices
        The puncture set ``S`` (column ``c`` <-> point ``c-1`` in ``F_p^m``).
    jmax : int or None
        Largest affine-span dimension to enumerate.  ``None`` -> ``m`` (EXACT).  A smaller
        cap returns the contribution of codewords with affine span ``<= jmax`` (a certified
        lower bound on ``A_d``; exact iff no higher-span codeword contributes).
    distance : int or None
        Target weight ``d``.  ``None`` -> computed via the certified MITM minimum distance.

    Returns
    -------
    dict with keys
        ``distance``            : the quantum distance d.
        ``A_d``                 : structured count of weight-d codewords (= true A_d when exact).
        ``low_weight_histogram``: {RM-weight w : count} of the contributing codewords
                                  (w = d + |supp(c) cap S|).
        ``dim_histogram``       : {affine-span dim j : count} of the contributing codewords.
        ``jmax``                : the cap actually used.
        ``exact``               : True iff jmax >= m (the full decomposition; no tail possible).
    """
    assert isinstance(p, int) and p >= 2
    rtilde = m * (p - 1) - r - 1
    code = build_triorthogonal_code(p, m, r, puncture_columns)
    if not code["full_rank"]:
        raise ValueError("puncture submatrix not full rank; the RM<->G0^perp bijection fails")
    G0 = np.asarray(code["G0"]).astype(np.int64) % p

    d = int(distance) if distance is not None else _distance(G0, p)
    if jmax is None:
        jmax = m
    jmax = min(jmax, m)

    S = set(int(c) - 1 for c in puncture_columns)
    pts = points(m, p)                                   # (p^m, m) for span computation
    A_d = 0
    rmweight_hist = Counter()
    dim_hist = Counter()
    found = [] if return_codewords else None

    for j in range(2, jmax + 1):
        rr = rtilde - (m - j) * (p - 1)
        if rr < 0:
            continue

        # --- fast path: j == 2, restricted code is the constants -> flat indicators ---
        if j == 2 and rr == 0:
            for colmap in _flats(m, j, p):
                surv = sum(1 for c in colmap if int(c) not in S)
                if surv == d:
                    A_d += (p - 1)
                    rmweight_hist[len(colmap)] += (p - 1)
                    dim_hist[2] += (p - 1)
                    if return_codewords:
                        for lam in range(1, p):
                            found.append({int(c): lam for c in colmap})
            continue

        Gj = np.asarray(rm_generator(rr, j, p)).astype(np.int64) % p   # dimCF x p^j
        dimCF = Gj.shape[0]
        # sound pruning: a span-j contributing codeword has |supp cap S| = w - d >=
        # min_span_weight(j) - d, so flats with fewer puncture points cannot host one.
        min_inter = max(0, _min_span_weight(rr, j, p) - d)

        for colmap in _flats(m, j, p):
            in_S = np.fromiter((int(c) in S for c in colmap), dtype=bool, count=len(colmap))
            n_inter = int(in_S.sum())
            if n_inter < min_inter:
                continue
            survpos = np.nonzero(~in_S)[0]
            if survpos.size < d:
                continue
            Gsurv = Gj[:, survpos]                         # dimCF x |survF|
            H = _null_space(Gsurv, p)                      # parity check of C_F|survF
            if H.shape[0] == 0:
                # survivors exactly span C_F: restricted code is the full space F_p^{|survF|}
                cws = _full_space_weight_d(int(survpos.size), d, p)
            else:
                cws = _mitm(H, d, p)
            if not cws:
                continue
            # one systematic decoder per flat (replaces the per-codeword linear solve)
            pivots, T = _systematic_lift(Gsurv, Gj, p)     # full = c[pivots] @ T
            piv_row = {int(c): r for r, c in enumerate(pivots)}
            for (sup_local, vals) in cws:
                # sup_local: indices into survpos; the codeword is determined by its values on
                # the info columns (pivots), so gather those and lift via T.
                cp = np.zeros(dimCF, dtype=np.int64)
                for i, v in zip(sup_local, vals):
                    r = piv_row.get(int(i))
                    if r is not None:
                        cp[r] = v % p
                full = (cp @ T) % p
                supp = np.nonzero(full)[0]                  # positions in F_p^j
                glob = colmap[supp]                         # global point indices
                w = int(supp.size)
                gp = pts[glob]
                diffs = (gp[1:] - gp[0]) % p
                adim = _rank(diffs, p) if diffs.shape[0] else 0
                if adim != j:
                    continue                                # minimal-flat dedup
                A_d += 1
                rmweight_hist[w] += 1
                dim_hist[j] += 1
                if return_codewords:
                    found.append({int(g): int(full[s]) for g, s in zip(glob, supp)})

    result = {
        "distance": d,
        "A_d": A_d,
        "low_weight_histogram": dict(sorted(rmweight_hist.items())),
        "dim_histogram": dict(sorted(dim_hist.items())),
        "jmax": jmax,
        "exact": jmax >= m,
    }
    if return_codewords:
        result["codewords"] = found
    return result
