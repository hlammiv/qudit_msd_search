"""Exact single-round distillation performance: leading order, acceptance, threshold.

Noise model (arXiv:2610.08012 sec 8): each input is (1-eps)|M><M| + eps/(p-1) sum_{j!=0}
Z^j|M><M|Z^-j, i.e. twirled depolarizing noise with the repo's ``delta`` = p/(p-1) eps.
Post-selecting on the X-stabilizers A = rowspan(G0) accepts Z(j) iff j in A^perp; output i
is wrong iff in addition j . h_i != 0, with h_i the i-th logical row.

Everything is driven by PRIMAL weight enumerators of small spaces C (dim C = dim G0 + 0/1/k),
using the per-coordinate MacWilliams collapse

    P(j in C^perp) = |C|^-1 sum_x C_x (1 - eta)^x ,   eta = p eps / (p-1).

Logical probabilities are differences of two such sums; we combine them into ONE integer
polynomial before evaluating (in mpmath at adaptive precision), so the ~eps^d logical part is
never formed by subtracting two O(1) floats.  Low-weight dual counts (d and the leading
coefficient) come from the same integer polynomial via Krawtchouk transforms at those
weights only -- the full O(n^3) MacWilliams transform is never needed.

Two output-error notions (identical for k = 1):
  * ``block``: P(any logical error | accept)      -- needs the enumerator of H = A + <H1>.
  * ``avg``:   (1/k) sum_i P(error on output i | accept) -- needs A + <h_i> for each i.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction

import mpmath
import numpy as np

from .weightdist import coset_weight_enumerators, krawtchouk, weight_enumerator


# ---------------------------------------------------------------------------
# Conventions
# ---------------------------------------------------------------------------
def eps_from_delta(delta, p):
    """Paper eps (prob. of a nontrivial Z^j) from the repo's depolarizing delta."""
    return (p - 1) / p * delta


def delta_from_eps(eps, p):
    return p / (p - 1) * eps


# ---------------------------------------------------------------------------
# Integer-polynomial representation
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class PrimalPoly:
    """P(eps) = scale * sum_x D_x (1 - eta)^x with integer D and rational scale.

    The dual (error-weight) count at weight w is scale * p^?-free Krawtchouk transform:
    N_w = scale * sum_x D_x K_w(x; n, p) -- an exact non-negative integer for every
    quantity built here (accepted errors, logical errors).
    """

    n: int
    p: int
    D: tuple[int, ...]
    scale: Fraction

    def __call__(self, eps, dps: int | None = None):
        p = self.p
        eta = mpmath.mpf(p) * mpmath.mpf(eps) / (p - 1)
        x = 1 - eta
        acc = mpmath.mpf(0)
        xp = mpmath.mpf(1)
        for Dx in self.D:
            if Dx:
                acc += Dx * xp
            xp *= x
        return acc * self.scale.numerator / self.scale.denominator

    def dual_count(self, w: int) -> int:
        """Number of error vectors of weight w counted by this polynomial (exact)."""
        tot = sum(Dx * krawtchouk(w, x, self.n, self.p) for x, Dx in enumerate(self.D) if Dx)
        val = self.scale * tot
        assert val.denominator == 1 and val >= 0, f"dual count at w={w} not a non-negative int: {val}"
        return int(val)


def _poly_sub(a, b):
    L = max(len(a), len(b))
    a = list(a) + [0] * (L - len(a))
    b = list(b) + [0] * (L - len(b))
    return tuple(x - y for x, y in zip(a, b))


# ---------------------------------------------------------------------------
# Profile
# ---------------------------------------------------------------------------
@dataclass
class DistillationProfile:
    """Exact performance data of one code: acceptance + logical polynomials."""

    p: int
    n: int
    k: int
    accept: PrimalPoly
    logical: dict = field(default_factory=dict)  # kind ("block"/"avg") -> PrimalPoly
    rigorous: bool = True
    method: str = "exact primal enumeration"

    def p_accept(self, eps):
        return self.accept(eps)

    def eps_out(self, eps, kind="avg"):
        with mpmath.workdps(self._dps(eps)):
            return self.logical[kind](eps) / self.accept(eps)

    def _dps(self, eps):
        # cancellation in sum_x D_x x^x costs ~log10(sum|D|) + d*log10(1/eps) digits
        mag = max(len(str(sum(abs(v) for v in P.D))) for P in self.logical.values())
        lost = 0 if eps >= 1 else int(math.ceil(self.n * 0 + 12 * max(1.0, -math.log10(max(eps, 1e-300)))))
        return 30 + mag + lost


def canonical_logical_rows(G, cols0, p):
    """Logical rows h_i = (u_i G) restricted to the survivors, with (u_i G)[cols0[j]] = delta_ij.

    For a triorthogonal RM space punctured at ``cols0`` this is the basis that makes the
    cubic form diagonal on the logicals (Lemma 11 generalised to k punctures), so it is the
    physical output basis.  Unique modulo rowspan(G0).
    """
    from .field import GFp

    GF = GFp(p)
    G = GF(np.asarray(G) % p)
    cols0 = list(cols0)
    sub = G[:, cols0]                       # (K, k)
    rref = sub.T.row_reduce()               # (k, K): pivots pick k independent rows of sub
    piv = []
    for i in range(rref.shape[0]):
        nz = np.flatnonzero(np.asarray(rref[i]))
        if nz.size:
            piv.append(int(nz[0]))
    assert len(piv) == len(cols0), "punctured-column submatrix not full rank"
    U = GF.Zeros((len(cols0), G.shape[0]))
    U[:, piv] = np.linalg.inv(sub[piv, :])
    full = U @ G
    keep = np.ones(G.shape[1], dtype=bool)
    keep[cols0] = False
    assert np.array_equal(np.asarray(full[:, cols0]), np.eye(len(cols0), dtype=np.int64))
    return np.asarray(full[:, keep], dtype=np.int64) % p


def profile_exact(G0, H1, p, kinds=("avg", "block"), max_words=50_000_000) -> DistillationProfile:
    """Exact profile from stabilizer rows ``G0`` (A) and logical rows ``H1``.

    ``avg`` costs p^dim(A) words (all k cosets A + h_i share one enumeration -- scalar
    multiples t h_i give the same weight histogram as h_i, so one coset per output suffices);
    ``block`` costs p^(dim A + k).  Raises ValueError if a requested kind exceeds max_words.
    """
    G0 = np.asarray(G0, dtype=np.int64) % p
    H1 = np.asarray(H1, dtype=np.int64) % p
    dA, k = G0.shape[0], H1.shape[0]

    offsets = [H1[i] for i in range(k)] if "avg" in kinds else []
    hists = coset_weight_enumerators(G0, [None] + offsets, p, max_words=max_words)
    Hw = None
    if "block" in kinds:
        Hw = weight_enumerator(np.vstack([G0, H1]), p, max_words=max_words)
    return profile_from_histograms(p, dA, hists[0], hists[1:] if "avg" in kinds else None, Hw)


def profile_from_histograms(p, dim_A, Aw, cosets=None, Hw=None) -> DistillationProfile:
    """Profile from primal histograms: ``Aw`` of A, ``cosets[i]`` of A + h_i, ``Hw`` of H.

    Engine-agnostic: the histograms may come from ``weightdist`` (local numpy) or from the
    ``dwq_coset`` C kernel.  ``cosets`` may cover only a SUBSET of the outputs (sampling);
    the "avg" polynomial is then the average over that subset, and k_avg records its size.
    """
    n = len(Aw) - 1
    sizeA = p ** dim_A
    assert sum(Aw) == sizeA, "A histogram does not sum to p^dim_A (rows dependent?)"
    accept = PrimalPoly(n, p, tuple(int(a) for a in Aw), Fraction(1, sizeA))
    logical = {}
    k = None
    if cosets is not None:
        k = len(cosets)
        assert all(sum(c) == sizeA for c in cosets), "coset histogram has wrong total"
        # P(err on i) = P(A^perp) - P((A+<h_i>)^perp) = (p-1)/(p|A|) sum_x (A_x - C^i_x) x^x
        Csum = [sum(int(c[x]) for c in cosets) for x in range(n + 1)]
        D = tuple(k * int(a) - c for a, c in zip(Aw, Csum))
        logical["avg"] = PrimalPoly(n, p, D, Fraction(p - 1, k * p * sizeA))
    if Hw is not None:
        kb = round(math.log(sum(Hw) // sizeA, p))
        assert sum(Hw) == sizeA * p ** kb, "rows of [G0; H1] must be independent"
        # P(A^perp) - P(H^perp) = |H|^-1 sum_x (p^k A_x - H_x) x^x
        D = _poly_sub([p ** kb * int(a) for a in Aw], [int(h) for h in Hw])
        logical["block"] = PrimalPoly(n, p, D, Fraction(1, sizeA * p ** kb))
        k = k if k is not None else kb
    return DistillationProfile(p=p, n=n, k=k or 0, accept=accept, logical=logical)


# ---------------------------------------------------------------------------
# Figures of merit
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class LeadingOrder:
    d: int
    L_d: Fraction          # weight-d logical error count (per output for "avg"; may be fractional)
    c_eps: float           # eps_out ~ c_eps * eps^d
    c_delta: float         # delta_out ~ c_delta * delta^d  (repo convention)


def leading_order(profile: DistillationProfile, kind="avg", wmax=None) -> LeadingOrder:
    """Smallest weight d with a logical accepted error, and the eps^d coefficient."""
    P = profile.logical[kind]
    p, n = profile.p, profile.n
    for w in range(1, (wmax or n) + 1):
        tot = sum(Dx * krawtchouk(w, x, n, p) for x, Dx in enumerate(P.D) if Dx)
        Lw = P.scale * tot
        if Lw:
            c_eps = float(Lw) / (p - 1) ** w
            c_delta = float(Lw) / ((p - 1) * p ** (w - 1))
            return LeadingOrder(d=w, L_d=Lw, c_eps=c_eps, c_delta=c_delta)
    raise ValueError("no logical error up to wmax")


@dataclass(frozen=True)
class Threshold:
    eps: float | None
    delta: float | None
    status: str            # "found" / "none below (p-1)/p"


def threshold(profile: DistillationProfile, kind="avg", grid=400, tol=1e-12) -> Threshold:
    """Smallest eps in (0, (p-1)/p) with eps_out(eps) = eps.

    eps = (p-1)/p (fully depolarized input) is always a trivial fixed point and is excluded.
    A log-spaced grid brackets the first sign change of eps_out - eps, then bisection.
    """
    p = profile.p
    top = (p - 1) / p
    f = lambda e: float(profile.eps_out(e, kind)) - e
    pts = list(np.logspace(-8, math.log10(top * (1 - 1e-6)), grid))
    prev_e, prev_f = pts[0], f(pts[0])
    for e in pts[1:]:
        fe = f(e)
        if prev_f < 0 <= fe:
            lo, hi = prev_e, e
            while hi - lo > tol * hi:
                mid = 0.5 * (lo + hi)
                if f(mid) < 0:
                    lo = mid
                else:
                    hi = mid
            root = 0.5 * (lo + hi)
            return Threshold(eps=root, delta=delta_from_eps(root, p), status="found")
        prev_e, prev_f = e, fe
    return Threshold(eps=None, delta=None, status="none below (p-1)/p")


# ---------------------------------------------------------------------------
# Catalog records -> matrices, tiers, evaluation
# ---------------------------------------------------------------------------
PERF_SCHEMA = "qmsd.performance.v1"
PERF_DIR = __import__("pathlib").Path(__file__).resolve().parent / "data" / "performance"

TIER1_WORDS = 3e8    # local numpy budget, in (codewords x cosets)
TIER2_WORDS = 1e13   # distributed C core (dist_weightdist) budget


MANHATTAN_EXPLICIT_MAX = 10 ** 4   # largest p^m for which a Manhattan record is rebuilt explicitly


def puncture_columns_of(record):
    """1-indexed puncture columns of a catalog record, or None if it has no explicit form.

    Stored columns win.  A small Manhattan record (``w`` set, p^m <= MANHATTAN_EXPLICIT_MAX) uses
    its weight-<=w point set.  Otherwise k <= 2 records use columns 1..k: AGL(m, p) is
    2-transitive on points, so any k <= 2 points are equivalent.
    """
    from .puncture import manhattan_puncture_columns

    if record.puncture_columns is not None:
        return tuple(record.puncture_columns)
    if record.m is None:
        return None
    if record.w is not None and record.p ** record.m <= MANHATTAN_EXPLICIT_MAX:
        cols = tuple(manhattan_puncture_columns(record.m, record.w, record.p))
        return cols if len(cols) == record.k else None
    if record.k <= 2:
        return tuple(range(1, record.k + 1))
    return None


def record_matrices(record):
    """(G0, H1) for a punctured-RM catalog record (columns from ``puncture_columns_of``)."""
    from .reedmuller import r_max, rm_generator
    from .triorthogonal import build_triorthogonal_code

    p, m = record.p, record.m
    assert m is not None, "record has no RM parameters"
    r = record.r if record.r is not None else r_max(m, p)
    cols = puncture_columns_of(record)
    assert cols is not None, "no explicit puncture set for this record"
    G = rm_generator(r, m, p)
    built = build_triorthogonal_code(p, m, r, list(cols), G=G)
    assert built["full_rank"], "puncture columns not full rank"
    G0 = np.asarray(built["G0"], dtype=np.int64) % p
    H1 = canonical_logical_rows(G, [c - 1 for c in cols], p)
    return G0, H1


def dim_stabilizer(record) -> int | None:
    from .reedmuller import monomial_exponents, r_max

    if record.m is None:
        return None
    r = record.r if record.r is not None else r_max(record.m, record.p)
    return len(monomial_exponents(record.m, r, record.p)) - record.k


def assign_tier(record) -> dict:
    """Which engine can handle each output-error kind for this record (see plan tiers 1-4)."""
    if puncture_columns_of(record) is None:
        return {"tier": 4, "dim_A": None, "engine_avg": "analytic", "engine_block": "analytic"}
    dA = dim_stabilizer(record)
    sizeA = record.p ** dA
    cost = {"avg": sizeA * max(1, record.k), "block": sizeA * record.p ** record.k}
    out = {"dim_A": dA, "log10_cost_avg": math.log10(cost["avg"]),
           "log10_cost_block": math.log10(cost["block"])}
    for kind, c in cost.items():
        out["engine_" + kind] = ("exact" if c <= TIER1_WORDS else
                                 "distributed" if c <= TIER2_WORDS else "bounds")
    out["tier"] = {"exact": 1, "distributed": 2, "bounds": 3}[out["engine_avg"]]
    return out


def summarize(profile: DistillationProfile, kind: str) -> dict:
    lo = leading_order(profile, kind)
    th = threshold(profile, kind)
    return {
        "d": lo.d,
        "L_d": str(lo.L_d),
        "c_eps": lo.c_eps,
        "c_delta": lo.c_delta,
        "threshold_eps": th.eps,
        "threshold_delta": th.delta,
        "threshold_status": th.status,
        "threshold_rigorous": profile.rigorous,
    }


def evaluate(record, tier_max: int = 1) -> dict:
    """Leading order + threshold for a catalog record, for every kind the tier allows."""
    t = assign_tier(record)
    res = {"artifact_id": record.artifact_id, "label": record.label, "p": record.p,
           "n": record.n, "k": record.k, "convention": "eps = paper (2610.08012); delta = p/(p-1) eps",
           **{k_: v for k_, v in t.items()}}
    kinds = tuple(k_ for k_ in ("avg", "block") if t.get("engine_" + k_) == "exact")
    if not kinds or t["tier"] > tier_max:
        res["status"] = "skipped (tier %s)" % t["tier"]
        return res
    G0, H1 = record_matrices(record)
    budget = int(TIER1_WORDS * max(1, record.p ** record.k))
    prof = profile_exact(G0, H1, record.p, kinds=kinds, max_words=budget)
    for kind in kinds:
        res[kind] = summarize(prof, kind)
    res["P_acc_eps_0.01"] = float(prof.p_accept(0.01))
    res["status"] = "exact"
    return res


def write_performance(res: dict, directory=None):
    """Store one ``evaluate`` result as ``<artifact_id>.json`` (schema qmsd.performance.v1)."""
    import json
    import subprocess
    from pathlib import Path

    d = Path(directory) if directory else PERF_DIR
    d.mkdir(parents=True, exist_ok=True)
    try:
        rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                             text=True, cwd=Path(__file__).parent).stdout.strip()
    except OSError:
        rev = ""
    out = {"schema": PERF_SCHEMA, "git_rev": rev, **res}
    path = d / f"{res['artifact_id']}.json"
    path.write_text(json.dumps(out, indent=1) + "\n")
    return path


def load_performance(artifact_id: str, directory=None) -> dict | None:
    import json
    from pathlib import Path

    path = (Path(directory) if directory else PERF_DIR) / f"{artifact_id}.json"
    return json.loads(path.read_text()) if path.exists() else None


# ---------------------------------------------------------------------------
# Tier 3: rigorous threshold lower bound from an enumerable subcode (paper Lemma 17)
# ---------------------------------------------------------------------------
def threshold_lower_bound_subcode(G0, H1, p, kind="avg", max_words=2e7, pool=24,
                                  ref_eps=None, seed=0) -> dict:
    """Rigorous lower bound on the threshold via a subcode A' of A = rowspan(G0).

    Lemma 17 (arXiv:2610.08012, FKG): with the logical rows fixed, A' <= A implies
    eps_out(A') >= eps_out(A) for all eps <= (p-1)/p, for any logical subspace -- so for both
    the block and each per-output error, hence their average.  Therefore
    threshold(A) >= threshold(A').  A' is grown greedily, one row at a time from a pool of
    random combinations of G0's rows, choosing the row that minimises eps_out(A') at
    ``ref_eps`` (default: the current bound), until p^dim(A') * cost would exceed max_words.
    """
    from .field import GFp

    GF = GFp(p)
    rng = np.random.default_rng(seed)
    G0 = np.asarray(G0, dtype=np.int64) % p
    H1 = np.asarray(H1, dtype=np.int64) % p
    k = H1.shape[0]
    per_word = k if kind == "avg" else p ** k
    rows: list[np.ndarray] = []
    best_thr = 0.0
    history = []
    while len(rows) < G0.shape[0] and p ** (len(rows) + 1) * per_word <= max_words:
        cands = (rng.integers(0, p, size=(pool, G0.shape[0])) @ G0) % p
        best = None
        for c in cands:
            trial = np.array(rows + [c])
            if np.linalg.matrix_rank(GF(np.vstack([trial, H1]))) < len(trial) + k:
                continue  # dependent over F_p (or hits the logical span)
            prof = profile_exact(trial, H1, p, kinds=(kind,), max_words=int(max_words))
            e = ref_eps or max(best_thr, 1e-3)
            val = float(prof.eps_out(e, kind))
            if best is None or val < best[0]:
                best = (val, c, prof)
        if best is None:
            break
        rows.append(best[1])
        th = threshold(best[2], kind)
        best_thr = max(best_thr, th.eps or 0.0)
        history.append({"dim": len(rows), "threshold_eps": th.eps})
    return {"threshold_eps_lower": best_thr or None,
            "threshold_delta_lower": delta_from_eps(best_thr, p) if best_thr else None,
            "subcode_dim": len(rows), "dim_A": int(G0.shape[0]), "history": history,
            "subcode": [r.tolist() for r in rows],
            "method": "Lemma 17 subcode bound (rigorous)"}
