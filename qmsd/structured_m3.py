"""Bounded m=3 structured-distance probe for the punctured-RM / triorthogonal codes.

Goal (see ``M3_PROBE_SPEC.md``): decide whether ``m=3`` can give *small* ``gamma<1``
codes at ``p=11/13`` -- specifically, whether a flat-spread puncture set can make the
**line + plane** structured distance clear ``gamma<1``, or whether the plane (2-flat)
codewords cap it (the ``m=3`` analog of the ``m=2`` 2D cap).

For a punctured triorthogonal code the quantum code is ``G0^perp`` = the *punctured*
generalized RM code ``RM_p(rtilde, m)`` restricted to the surviving coordinates
(``rtilde = m(p-1) - r - 1``).  When the punctured-column submatrix is full rank the
restriction map ``RM_p(rtilde,m) -> G0^perp`` is a bijection, so

    d = min_{c != 0 in RM_p(rtilde,m)} |supp(c) \\ S|        (S = puncture set, |S| = k).

The low-weight codewords are geometrically structured (Delsarte-Goethals-MacWilliams,
Kasami-Tokura, Leducq): every codeword is supported on an affine flat ``F`` and restricts
to a codeword of ``C_F = RM_p(rtilde - (m-j)(p-1), j)`` on ``F ~= F_p^j``.  This module
covers the ``j = 1`` (line) and ``j = 2`` (plane) classes for ``m = 3``:

    d_lines  = d_RM - max_{line ell}  |S cap ell|     (1-flat; closed form, beta_line=0)
    d_planes = min over planes F of the min punctured weight of a plane-supported codeword
    d_struct = min(d_lines, d_planes)                 (UPPER bound on the true distance d)

``d_struct >= d_true`` because it minimises over a *subset* of codewords (flats of dim<=2);
the dim-3 full-span codewords -- deferred to the full build -- can only lower it further.
So if even ``d_struct`` cannot clear ``gamma<1`` the answer is a clean NO-GO.

For ``p=11, m=3``: ``r_max=9, rtilde=20, d_RM=11``; the plane code is ``RM_11(10,2)`` with
``beta=0`` (genuinely-2D minimum weight ``w2 = 2(p-1) = 20``).  ``p=13`` is analogous
(``rtilde=24``, plane code ``RM_13(12,2)``, ``w2=24``).

Point-index convention matches ``qmsd.reedmuller.points``: column ``c`` <-> point
``(x,y,z) = (c%p, (c//p)%p, (c//p^2)%p)``, i.e. ``c = x + p*y + p^2*z`` (x least significant).
"""
from __future__ import annotations

import itertools

import numpy as np

from .reedmuller import r_tilde, d_rm, points, rm_generator
from .structured_ad import _rank, _null_space
from .triorthogonal import build_triorthogonal_code

__all__ = [
    "enumerate_lines_3d",
    "enumerate_planes_3d",
    "line_distance",
    "plane_distance",
    "flat_distance",
    "flat_spread_puncture",
]


# ---------------------------------------------------------------------------
# affine flats of AG(3,p)
# ---------------------------------------------------------------------------
def _canonical_directions(p):
    """Yield each 1-dim direction of F_p^3 once, as a length-3 int tuple with leading 1."""
    # leading nonzero coordinate = index of first 1; entries before it are 0.
    for a in range(p):
        for b in range(p):
            yield (1, a, b)
    for b in range(p):
        yield (0, 1, b)
    yield (0, 0, 1)


def enumerate_lines_3d(p):
    """Yield ``(point_indices, t_params)`` for every affine line of AG(3,p).

    A line is ``{a + t*v : t in F_p}`` for a canonical direction ``v`` (leading nonzero
    coordinate = 1) and a base point ``a`` taken from the transversal hyperplane
    ``a[lead] = 0``.  Each of the ``(p^2+p+1)*p^2`` lines is produced exactly once.
    ``point_indices`` are the ``p`` global column indices; ``t_params`` is ``[0..p-1]`` so a
    univariate codeword in ``t`` lifts directly (matches ``structured_pe`` convention).
    """
    POW = np.array([1, p, p * p], dtype=np.int64)
    Fp = range(p)
    rng = list(Fp)
    for v in _canonical_directions(p):
        v = np.array(v, dtype=np.int64)
        lead = int(np.nonzero(v)[0][0])
        # base points a with a[lead] = 0
        free = [i for i in range(3) if i != lead]
        for vals in itertools.product(Fp, repeat=2):
            a = np.zeros(3, dtype=np.int64)
            a[free[0]] = vals[0]
            a[free[1]] = vals[1]
            pts = (a[None, :] + np.outer(np.array(rng), v)) % p  # (p,3)
            idx = (pts @ POW).astype(np.int64)
            yield idx.tolist(), rng


def _canonical_normals(p):
    """Yield each hyperplane normal direction of F_p^3 once (leading nonzero = 1)."""
    yield from _canonical_directions(p)


def enumerate_planes_3d(p):
    """Yield ``(point_indices, colmap_inv)`` for every affine plane (2-flat) of AG(3,p).

    A plane is the affine hyperplane ``{x : n . x = c}`` for a canonical normal ``n`` and
    ``c in F_p``; there are ``(p^2+p+1)*p = 1463`` of them for ``p=11``.  ``point_indices`` is
    the array of the ``p^2`` global columns of the plane, ordered so that internal AG(2,p)
    coordinate ``(u,v)`` (index ``u + p*v``) maps to ``point_indices[u + p*v]`` -- i.e. the
    columns align with ``rm_generator(*, 2, p)`` and ``points(2,p)``.
    """
    POW = np.array([1, p, p * p], dtype=np.int64)
    Fp = list(range(p))
    pts2 = points(2, p)  # (p^2, 2): internal coords (u,v), u least significant
    for n in _canonical_normals(p):
        n = np.array(n, dtype=np.int64)
        lead = int(np.nonzero(n)[0][0])
        # two basis vectors b1,b2 spanning the linear hyperplane n.x=0
        basis = []
        for i in range(3):
            if i == lead:
                continue
            b = np.zeros(3, dtype=np.int64)
            b[i] = 1
            # set lead coord so that n.b = 0: n[lead]*b[lead] + n[i]*1 = 0
            inv = pow(int(n[lead]), p - 2, p)
            b[lead] = (-n[i] * inv) % p
            basis.append(b)
        b1, b2 = basis
        for c in Fp:
            # an origin o with n.o = c: put c/n[lead] on the lead coord
            o = np.zeros(3, dtype=np.int64)
            o[lead] = (c * pow(int(n[lead]), p - 2, p)) % p
            # point(u,v) = o + u*b1 + v*b2
            coords = (o[None, :] + pts2[:, 0:1] * b1[None, :] + pts2[:, 1:2] * b2[None, :]) % p
            idx = (coords @ POW).astype(np.int64)
            yield idx, np.arange(p * p)  # colmap is identity: internal index == row index


# ---------------------------------------------------------------------------
# line-supported distance  (1-flats of AG(3,p))
# ---------------------------------------------------------------------------
def line_distance(p, m, r, puncture_columns_1indexed, G=None, X=None, n=None):
    """Line-supported punctured distance ``d_lines`` of the punctured RM code (m=3).

    On a line ``ell`` the codewords of ``RM_p(rtilde,m)`` supported in ``ell`` are the
    univariate polynomials of degree ``<= beta_line := rtilde - (m-1)(p-1)`` (the restricted
    code ``RM_p(beta_line, 1)``).  A degree-``<=beta_line`` poly has ``<= beta_line`` roots, so
    on a line with ``s = |S cap ell|`` punctures the minimum surviving weight is
    ``(p - s) - beta_line``; hence

        d_lines = (p - beta_line) - max_ell |S cap ell| = d_RM - max_ell |S cap ell|.

    Returns ``{d_lines, d_RM, beta_line, max_line_punctures, best_line_points, certificate}``.
    The certificate is an explicit weight-``d_lines`` codeword verified in ``ker(X_stab)``.
    """
    assert m == 3, "line_distance here is specialised to m=3 (AG(3,p))"
    rtilde = r_tilde(r, m, p)
    beta_line = rtilde - (m - 1) * (p - 1)
    if beta_line < 0:
        raise NotImplementedError(f"beta_line={beta_line} < 0; line-supported structure n/a")
    d_RM = d_rm(rtilde, m, p)
    # for the p=11/13 m=3 codes beta_line == 0 and d_RM == p; assert the identity generally
    assert d_rm(beta_line, 1, p) == p - beta_line

    S = set(int(c) - 1 for c in puncture_columns_1indexed)

    max_pun = -1
    best_line = None
    for line_pts, t_params in enumerate_lines_3d(p):
        s = sum(1 for pt in line_pts if pt in S)
        if s > max_pun:
            max_pun = s
            best_line = (line_pts, t_params)

    d_lines = (p - beta_line) - max_pun
    assert d_lines == d_RM - max_pun

    cert = None
    if X is not None and n is not None:
        cert = _line_certificate(best_line, S, beta_line, p, n, X, puncture_columns_1indexed)

    return {
        "d_lines": int(d_lines),
        "d_RM": int(d_RM),
        "beta_line": int(beta_line),
        "max_line_punctures": int(max_pun),
        "best_line_points": list(best_line[0]),
        "certificate": cert,
    }


# ---------------------------------------------------------------------------
# plane-supported distance  (2-flats of AG(3,p))
# ---------------------------------------------------------------------------
def plane_distance(p, m, r, puncture_columns_1indexed, d_max=6, G0=None, S=None,
                   only_heavy=None, mdc=None, return_witness=False):
    """Plane-supported punctured distance ``d_planes`` of the punctured RM code (m=3).

    For each affine plane (2-flat) ``F`` of AG(3,p), the codewords of ``G0^perp`` whose support
    lies inside ``F`` are exactly the (extended-by-zero) codewords of the punctured plane code
    ``Punc_{S cap F}( RM_p(rtilde-(m-2)(p-1), 2) )``; their punctured weight ``|supp(c)\\S|`` is
    the number of *surviving* columns of ``F`` they use.  That minimum is precisely the minimum
    number of F_p-dependent columns of ``G0`` restricted to ``F``'s surviving columns:

        d_planes = min over planes F of  min_dependent_columns( G0[:, surv(F)], p ).

    This is CERTIFIED (the meet-in-the-middle column-dependency search) and captures BOTH the
    1-flat (line) and 2-flat (plane) codewords -- a line lies inside a plane -- so the overall
    ``min(d_lines, d_planes)`` equals ``min over planes F of min_dependent_columns(G0[:,surv F])``.
    Only the genuinely 3-dimensional (full-span) codewords are excluded; they can only lower the
    true distance further, so this is an UPPER bound on ``d``.

    ``d_max`` caps the per-plane search depth (the MITM certifies up to 6).  ``only_heavy`` (an int)
    restricts the scan to planes with ``|S cap F| >= only_heavy`` (a sound prune only when the
    minimal plane codeword weight ``w2`` satisfies ``w2 - only_heavy >= d_max``; callers must
    justify it).  ``mdc`` overrides the min-dependent-columns routine (e.g. an overflow-safe one).
    """
    assert m == 3
    if mdc is None:
        from .mindist import min_dependent_columns as mdc
    if G0 is None:
        code = build_triorthogonal_code(p, m, r, puncture_columns_1indexed)
        if not code["full_rank"]:
            raise ValueError("puncture submatrix not full rank; RM<->G0^perp bijection fails")
        G0 = np.asarray(code["G0"]).astype(np.int64) % p
    else:
        G0 = np.asarray(G0).astype(np.int64) % p
    if S is None:
        S = set(int(c) - 1 for c in puncture_columns_1indexed)
    Sset = set(int(c) for c in S)

    surv = [c for c in range(p ** m) if c not in Sset]
    surv_index = {c: i for i, c in enumerate(surv)}

    best = d_max + 1
    witness = None
    n_scanned = 0
    for idxF, _cm in enumerate_planes_3d(p):
        idxF = idxF.tolist()
        nS = sum(1 for c in idxF if c in Sset)
        if only_heavy is not None and nS < only_heavy:
            continue
        cols = [surv_index[c] for c in idxF if c not in Sset]
        if len(cols) < 1:
            continue
        H = G0[:, cols]
        n_scanned += 1
        try:
            dF = mdc(H, p, d_max=min(d_max, best - 1) if best <= d_max else d_max)
        except ValueError:
            continue  # no dependent set <= cap on this plane
        if dF < best:
            best = dF
            if return_witness:
                witness = {"plane_points": idxF, "n_S_in_plane": nS, "d": int(dF)}
            if best <= 1:
                break
    return {
        "d_planes": int(best) if best <= d_max else None,
        "d_planes_capped_at": int(d_max),
        "exceeds_d_max": best > d_max,
        "planes_scanned": n_scanned,
        "witness": witness,
    }


def flat_distance(p, m, r, puncture_columns_1indexed, d_max=6, mdc=None):
    """Structured flat distance ``{d_lines, d_planes, d_struct, witness}`` of the m=3 code.

    ``d_struct = min(d_lines, d_planes)`` is an UPPER bound on the true distance ``d`` (it omits
    only the dim-3 full-span codewords).  ``d_planes`` is searched up to ``d_max`` (MITM cap 6);
    if every plane exceeds ``d_max`` it is reported as ``None`` with ``exceeds_d_max=True`` and
    ``d_struct`` falls back to ``min(d_lines, d_max+...)`` semantics handled by the caller.
    """
    assert m == 3
    code = build_triorthogonal_code(p, m, r, puncture_columns_1indexed)
    if not code["full_rank"]:
        raise ValueError("puncture submatrix not full rank; RM<->G0^perp bijection fails")
    G0 = np.asarray(code["G0"]).astype(np.int64) % p
    X = np.asarray(code["X_stab"]).astype(np.int64) % p
    n = code["n"]
    S = set(int(c) - 1 for c in puncture_columns_1indexed)

    lin = line_distance(p, m, r, puncture_columns_1indexed, X=X, n=n)
    pl = plane_distance(p, m, r, puncture_columns_1indexed, d_max=d_max, G0=G0, S=S,
                        mdc=mdc, return_witness=True)

    d_lines = lin["d_lines"]
    d_planes = pl["d_planes"]
    if d_planes is None:
        # every plane's min weight exceeds d_max; planes do not undercut d_max
        d_struct = min(d_lines, d_max + 1)
        d_struct_note = f"d_planes > {d_max} (not undercutting); d_struct = min(d_lines, >{d_max})"
    else:
        d_struct = min(d_lines, d_planes)
        d_struct_note = None
    return {
        "d_lines": d_lines,
        "d_planes": d_planes,
        "d_struct": int(d_struct),
        "d_max": d_max,
        "line_info": lin,
        "plane_info": pl,
        "note": d_struct_note,
        "dim_G0": int(G0.shape[0]),
        "n": n,
        "k": code["k"],
    }


# ---------------------------------------------------------------------------
# flat-spread puncture sampler (minimise max_flat |S cap flat| over lines & planes)
# ---------------------------------------------------------------------------
def _line_incidence(p):
    """(lines list, point->line-id arrays) for AG(3,p).  Memoised per p on the function."""
    cache = _line_incidence._cache
    if p in cache:
        return cache[p]
    lines = [idx for idx, _ in enumerate_lines_3d(p)]
    N = p ** 3
    pt_lines = [[] for _ in range(N)]
    for li, idx in enumerate(lines):
        for c in idx:
            pt_lines[c].append(li)
    pt_lines = [np.array(x, dtype=np.int64) for x in pt_lines]
    cache[p] = (lines, pt_lines)
    return cache[p]


_line_incidence._cache = {}


def flat_spread_puncture(p, k, seed=0, ls_iters=4000, cand=400, return_maxline=False):
    """A flat-spread size-``k`` puncture set of AG(3,p): minimise ``max_line |S cap line|``.

    Greedy incremental construction (each point added is the one with the smallest current
    maximum incident line-count) followed by min-max local search (swap a point off a heaviest
    line for one that does not raise the maximum).  A small ``max_line`` directly maximises
    ``d_lines = d_RM - max_line`` and, because every minimal (weight-``2(p-1)``) plane codeword
    is supported on a union of two lines, also bounds the minimal-2D plane punctures by
    ``2 max_line`` -- so a line-spread set is the natural flat-spread candidate.

    Returns the 0-indexed column array (and ``max_line`` if ``return_maxline``).
    """
    lines, pt_lines = _line_incidence(p)
    L = len(lines)
    N = p ** 3
    rng = np.random.default_rng(seed)
    line_count = np.zeros(L, dtype=np.int32)
    chosen = np.zeros(N, dtype=bool)
    avail = np.ones(N, dtype=bool)
    nchosen = 0
    while nchosen < k:
        pool = np.nonzero(avail)[0]
        if pool.size > cand:
            pool = rng.choice(pool, size=cand, replace=False)
        best_c = -1
        best_score = None
        for c in pool:
            mm = int(line_count[pt_lines[c]].max())
            if best_score is None or mm < best_score:
                best_score = mm
                best_c = c
                if mm == 0:
                    break
        chosen[best_c] = True
        avail[best_c] = False
        nchosen += 1
        line_count[pt_lines[best_c]] += 1
    ml = int(line_count.max())
    for _ in range(ls_iters):
        heavy = np.nonzero(line_count >= ml)[0]
        if heavy.size == 0:
            break
        idx = lines[heavy[rng.integers(heavy.size)]]
        on = [c for c in idx if chosen[c]]
        if not on:
            continue
        rem = on[rng.integers(len(on))]
        for _ in range(15):
            add = int(rng.integers(N))
            if chosen[add]:
                continue
            line_count[pt_lines[rem]] -= 1
            line_count[pt_lines[add]] += 1
            if int(line_count.max()) <= ml:
                chosen[rem] = False
                chosen[add] = True
                break
            line_count[pt_lines[rem]] += 1
            line_count[pt_lines[add]] -= 1
        ml = int(line_count.max())
    cols = np.nonzero(chosen)[0]
    return (cols, ml) if return_maxline else cols


def _line_certificate(best_line, S, beta_line, p, n, X, puncture_columns_1indexed):
    """Explicit weight-d_lines codeword on the heaviest line, verified in ker(X_stab)."""
    line_pts, t_params = best_line
    surv = [c for c in range(p ** 3) if c not in S]
    surv_index = {c: i for i, c in enumerate(surv)}
    surv_on = [(pt, t) for pt, t in zip(line_pts, t_params) if pt not in S]
    roots = [t for _pt, t in surv_on[:beta_line]]

    def g(t):
        v = 1
        for rt in roots:
            v = (v * ((t - rt) % p)) % p
        return v

    v = np.zeros(n, dtype=np.int64)
    for pt, t in zip(line_pts, t_params):
        if pt not in S:
            v[surv_index[pt]] = g(t) % p
    weight = int(np.count_nonzero(v))
    in_kernel = bool(np.all((np.asarray(X).astype(np.int64) @ v) % p == 0))
    return {"weight": weight, "in_kernel": in_kernel}
