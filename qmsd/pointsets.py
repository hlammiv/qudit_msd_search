"""Point-set tools for the low-weight words of RM_p(r, m)^perp = RM_p(m(p-1)-1-r, m).

A vector u supported on a point set T of F_p^m is a word of RM(r, m)^perp iff the columns of the
RM(r, m) evaluation matrix at T are dependent with coefficients u; u has FULL support on T iff no
column of T is a coloop (every point lies in a circuit of the column matroid).  Over F_p (p >= 3)
a vector space is never the union of <= p proper subspaces, so "no coloop" is equivalent to the
existence of a single null vector with all entries non-zero.

Used by the weight-5/6 line lemma and the weight-7 lemma of WEIGHTED_SUPPORTS.md (the mixed-word
stages of the p = 7 certificates): batched F_p rank over many small matrices, affine lines in
parameter order, and an affine-invariant classification of small point sets (collinearity
profile, span, maximum coplanar subset, quadric membership).
"""
from __future__ import annotations

import itertools

import numpy as np

from .reedmuller import points
from .structured_ad import _null_space, _rank, _rref


def inverse_table(p):
    """inv[a] = a^{-1} mod p (inv[0] = 0)."""
    return np.array([0] + [pow(a, -1, p) for a in range(1, p)], dtype=np.int64)


def batched_rank_modp(A, p):
    """Ranks over F_p of a stack of matrices A (B, n, m), by batched Gaussian elimination."""
    A = (np.asarray(A, dtype=np.int64) % p).copy()
    B, n, m = A.shape
    inv = inverse_table(p)
    rank = np.zeros(B, dtype=np.int64)
    row = np.zeros(B, dtype=np.int64)
    ar = np.arange(n)[None, :]
    for c in range(m):
        mask = (A[:, :, c] != 0) & (ar >= row[:, None])
        has = mask.any(axis=1)
        if not has.any():
            continue
        idx = np.nonzero(has)[0]
        r0 = row[idx]
        r1 = np.argmax(mask[idx], axis=1)
        sub = A[idx]                                   # (I, n, m) copy
        ii = np.arange(len(idx))
        tmp = sub[ii, r0].copy()
        sub[ii, r0] = sub[ii, r1]
        sub[ii, r1] = tmp
        piv = sub[ii, r0, c]
        sub[ii, r0] = (sub[ii, r0] * inv[piv][:, None]) % p
        fac = sub[:, :, c].copy()
        fac[ii, r0] = 0
        sub = (sub - fac[:, :, None] * sub[ii, r0][:, None, :]) % p
        A[idx] = sub
        row[idx] += 1
        rank[idx] += 1
        if (row >= n).all():
            break
    return rank


def full_support_null_vector(M, p):
    """A null vector of M (rows x n) over F_p with every entry non-zero, or None.

    Enumerates the null space when its dimension is <= 4 (p^4 combinations), otherwise falls back
    to the coloop test (column i is a coloop iff it is zero on the whole null space) and a random
    search.
    """
    M = np.asarray(M, dtype=np.int64) % p
    ns = np.asarray(_null_space(M, p)) % p
    if ns.shape[0] == 0:
        return None
    if (ns != 0).any(axis=0).sum() < M.shape[1]:          # a coordinate that is zero on the null space
        return None
    if ns.shape[0] <= 4:
        for coef in itertools.product(range(p), repeat=ns.shape[0]):
            if any(coef):
                u = (np.asarray(coef) @ ns) % p
                if u.all():
                    return u
        return None
    rng = np.random.default_rng(0)
    for _ in range(100000):
        u = (rng.integers(0, p, ns.shape[0]) @ ns) % p
        if u.all():
            return u
    raise RuntimeError("no coloop but no full-support vector found (should not happen for p >= 3)")


def affine_lines(p, m):
    """All affine lines of F_p^m, as an (L, p) array of point indices (``reedmuller.points``
    order), row l = a + t v for t = 0 .. p-1.  L = p^(m-1) (p^m - 1) / (p - 1)."""
    POW = p ** np.arange(m)
    dirs = [v for v in itertools.product(range(p), repeat=m)
            if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1]
    rows = []
    for v in dirs:
        piv = next(i for i, x in enumerate(v) if x)
        reps = points(m - 1, p)                                        # free coordinates except the pivot
        a = np.zeros((reps.shape[0], m), dtype=np.int64)
        a[:, [i for i in range(m) if i != piv]] = reps
        pts = (a[:, None, :] + np.arange(p)[None, :, None] * np.asarray(v)[None, None, :]) % p
        rows.append(pts @ POW)
    return np.concatenate(rows, axis=0)


def _affine_rank(pts, p):
    pts = np.asarray(pts, dtype=np.int64) % p
    if len(pts) <= 1:
        return 0
    return _rank((pts[1:] - pts[0]) % p, p)


def line_through(a, b, p):
    """The p points of the line through distinct a, b (as a sorted tuple of coordinate tuples)."""
    a = np.asarray(a) % p
    v = (np.asarray(b) - a) % p
    return tuple(sorted(tuple(int(x) for x in (a + t * v) % p) for t in range(p)))


def _solve_coords(Bm, d, p):
    """c with c @ Bm = d over F_p (Bm has independent rows)."""
    k = Bm.shape[0]
    aug = np.concatenate([Bm.T % p, (d % p)[:, None]], axis=1)       # (m, k+1)
    R, piv = _rref(aug, p)
    assert k not in piv, "inconsistent system"
    c = np.zeros(k, dtype=np.int64)
    for r, pc in enumerate(piv):
        c[pc] = R[r, k]
    return c


def collinearity_profile(pts, p):
    """Affine-invariant description of a small point set (coordinates, (n, m)):
    span, sizes of the lines with >= 3 points (descending), max collinear, max coplanar
    (for span >= 3), and whether the set lies on a quadric of its own span."""
    pts = np.asarray(pts, dtype=np.int64) % p
    n, m = pts.shape
    span = _affine_rank(pts, p)
    lines = {}
    for i, j in itertools.combinations(range(n), 2):
        L = line_through(pts[i], pts[j], p)
        lines.setdefault(L, set()).update((i, j))
    sizes = sorted((len(s) for s in lines.values()), reverse=True)
    max_col = sizes[0] if sizes else min(n, 1)
    profile = tuple(s for s in sizes if s >= 3)
    max_cop = None
    if span >= 3:
        max_cop = 0
        for tri in itertools.combinations(range(n), 3):
            if _affine_rank(pts[list(tri)], p) < 2:
                continue
            base = pts[tri[0]]
            B = (pts[[tri[1], tri[2]]] - base) % p
            cnt = sum(1 for x in range(n) if _rank(np.vstack([B, (pts[x] - base) % p]), p) <= 2)
            max_cop = max(max_cop, cnt)
    on_quadric = None                       # only meaningful when the set has >= C(span+2, 2) points
    if span >= 2 and n >= (span + 1) * (span + 2) // 2:
        base = pts[0]
        diffs = (pts - base) % p
        basis = []
        for d in diffs:
            if d.any() and _rank(np.vstack(basis + [d]), p) == len(basis) + 1:
                basis.append(d)
        Bm = np.vstack(basis)                                          # (span, m)
        coords = np.array([_solve_coords(Bm, d, p) for d in diffs])    # (n, span)
        exps = [e for e in itertools.product(range(3), repeat=span) if sum(e) <= 2]
        mono = np.array([[int(np.prod([pow(int(c), a, p) for c, a in zip(row, e)])) % p for e in exps]
                         for row in coords])
        on_quadric = bool(_rank(mono, p) < len(exps))
    return {"n": n, "span": span, "lines_ge3": profile, "max_collinear": max_col,
            "max_coplanar": max_cop, "on_quadric": on_quadric}


def classify_pointset(pts, p):
    """Short label: 'line7' for seven collinear points, 'line6+1', 'line5+2' (>= 5 collinear),
    otherwise 'span{j}:lines{profile}[:quadric]'."""
    prof = collinearity_profile(pts, p)
    n = prof["n"]
    if prof["max_collinear"] == n:
        return f"line{n}"
    if prof["max_collinear"] >= 5:
        return f"line{prof['max_collinear']}+{n - prof['max_collinear']}"
    lab = f"span{prof['span']}:lines{list(prof['lines_ge3'])}"
    if prof["on_quadric"]:
        lab += ":quadric"
    return lab


def frame_and_flat(m, j, p):
    """(frame, flat): point indices of 0, e_1..e_j and of all points of the coordinate j-flat
    {x : x_{j+1} = .. = x_m = 0}."""
    P = points(m, p)
    index = {tuple(int(v) for v in x): i for i, x in enumerate(P)}
    frame = [index[(0,) * m]] + [index[tuple(1 if i == a else 0 for i in range(m))] for a in range(j)]
    flat = [i for i, x in enumerate(P) if not x[j:].any()]
    return frame, flat


# ---------------------------------------------------------------------------
# direction tables for the "extend by <= 2 class columns" test of the mixed stages
# ---------------------------------------------------------------------------
def _chunk(p):
    """Coordinates per int64 key: p^c < 2^62."""
    c = 1
    while p ** (c + 1) < 2 ** 62:
        c += 1
    return c


def direction_keys(V, p):
    """Exact keys of the projective directions of the rows of V (zero rows -> all-zero key): the
    rows are scaled so that the first non-zero entry is 1 and read, in chunks of coordinates, as
    base-p integers.  Returns an (n, n_chunks) int64 array (one chunk for <= 22 coordinates at
    p = 7); equal rows <=> equal keys."""
    V = np.asarray(V, dtype=np.int64) % p
    inv = inverse_table(p)
    first = np.argmax(V != 0, axis=1)
    lead = V[np.arange(len(V)), first]
    Vn = (V * inv[lead][:, None]) % p
    c = _chunk(p)
    cols = []
    for lo in range(0, V.shape[1], c):
        blk = Vn[:, lo:lo + c]
        cols.append((blk @ (p ** np.arange(blk.shape[1])).astype(np.int64)).astype(np.int64))
    return np.stack(cols, axis=1)


def _lexsort_rows(K):
    return np.lexsort(tuple(K[:, i] for i in range(K.shape[1] - 1, -1, -1)))


def _find_rows(sorted_K, Q):
    """Positions in sorted_K (lexicographically sorted rows) of the rows of Q, or -1."""
    n = len(sorted_K)
    lo = np.zeros(len(Q), dtype=np.int64); hi = np.full(len(Q), n, dtype=np.int64)
    # refine the [lo, hi) range chunk by chunk: binary search within the previous range
    # (vectorised for the first chunk; the following chunks only act on the few multi-entry ranges)
    k0 = sorted_K[:, 0]
    lo = np.searchsorted(k0, Q[:, 0], "left"); hi = np.searchsorted(k0, Q[:, 0], "right")
    out = np.full(len(Q), -1, dtype=np.int64)
    one = hi - lo == 1
    out[one] = np.where((sorted_K[lo[one]] == Q[one]).all(axis=1), lo[one], -1)
    for qi in np.nonzero(hi - lo > 1)[0]:
        block = sorted_K[lo[qi]:hi[qi]]
        m = np.nonzero((block == Q[qi]).all(axis=1))[0]
        if len(m):
            out[qi] = lo[qi] + m[0]
    return out


def pair_direction_table(C, p):
    """Lookup table of every direction spanned by <= 2 columns of C (d x N): the zero vector, the
    N single columns and the (p-1) C(N,2) combinations c_i + beta c_j.  Returns a dict with the
    sorted int64 keys ('keys') and, aligned with them, 'kind' (0 zero, 1 single, 2 pair), the
    column indices 'i', 'j' (-1 if unused) and 'beta'.  A vector s is a combination of <= 2
    columns iff direction_keys(s) is in 'keys' (membership by np.searchsorted)."""
    C = np.asarray(C, dtype=np.int64) % p
    N = C.shape[1]
    single = C.T.copy()
    ii, jj = np.triu_indices(N, 1)
    pair_vecs = np.concatenate([(single[ii] + beta * single[jj]) % p for beta in range(1, p)], axis=0)
    keys = np.concatenate([direction_keys(np.zeros((1, C.shape[0]), dtype=np.int64), p),
                           direction_keys(single, p), direction_keys(pair_vecs, p)], axis=0)
    kind = np.concatenate([[0], np.ones(N, dtype=np.int64), np.full(len(pair_vecs), 2, dtype=np.int64)])
    i = np.concatenate([[-1], np.arange(N), np.tile(ii, p - 1)])
    j = np.concatenate([[-1], np.full(N, -1), np.tile(jj, p - 1)])
    beta = np.concatenate([[0], np.zeros(N, dtype=np.int64), np.repeat(np.arange(1, p), len(ii))])
    order = _lexsort_rows(keys)
    return {"keys": keys[order], "kind": kind[order], "i": i[order], "j": j[order], "beta": beta[order]}


def lookup_directions(table, S, p):
    """Positions in table['keys'] of the directions of the rows of S (d-vectors), or -1."""
    return _find_rows(table["keys"], direction_keys(S, p))
