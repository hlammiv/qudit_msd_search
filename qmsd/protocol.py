"""Multi-round distillation protocols from exact per-output curves (cost to reach a target error).

A protocol is a sequence of codes Q_1..Q_r.  Round i consumes n_i states at per-output error
eps_{i-1} and, with probability P_acc,i(eps_{i-1}), emits k_i states at per-output error
eps_i = f_i(eps_{i-1}) (``qmsd.performance`` "avg" curve).  Inputs of each round are drawn
from distinct blocks of the previous round, so the per-output marginal is the right input
model, and the output keeps the twirled form (arXiv:2610.08012 sec 8).  The expected cost in
raw magic states per final output is

    C = prod_i n_i / (k_i * P_acc,i(eps_{i-1})).

``best_protocol`` searches all sequences (with repetition) from a code library by
depth-first search with cost pruning.  Codes whose logical gate is U^lambda with lambda
not a cube (p = 1 mod 3) are excluded: their output is in a different Clifford class.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import mpmath
import numpy as np

from .performance import (PERF_DIR, DistillationProfile, leading_order, profile_from_histograms,
                          threshold)

HIST_DIR = PERF_DIR / "hist"
LOCAL_NUMPY_WORDS = 3e8   # above this, histograms come from the local dwq_coset C kernel
LOCAL_C_THREADS = 4
EXACT_BELOW = 1e-7        # below this eps, use the leading order (relative error ~ n*eps)


@dataclass
class CodeCurve:
    """Exact per-output single-round map of one code."""

    name: str
    p: int
    n: int
    k: int
    profile: DistillationProfile = field(repr=False)
    d: int = 0
    c_eps: float = 0.0
    threshold_eps: float | None = None
    note: str = ""

    def __post_init__(self):
        lo = leading_order(self.profile, "avg")
        self.d, self.c_eps = lo.d, lo.c_eps
        self.threshold_eps = threshold(self.profile, "avg").eps

    def f(self, eps):
        """Per-output output error (mpf; exact above EXACT_BELOW, leading order below)."""
        if eps < EXACT_BELOW:
            return mpmath.mpf(self.c_eps) * mpmath.mpf(eps) ** self.d
        return self.profile.eps_out(eps, "avg")

    def pacc(self, eps) -> float:
        if eps < EXACT_BELOW:
            return max(0.0, 1.0 - self.n * float(eps))
        with mpmath.workdps(30):
            return float(self.profile.p_accept(eps))

    @property
    def gamma(self) -> float:
        return math.log(self.n / self.k) / math.log(self.d) if self.d > 1 else float("inf")

    # bounded-curve interface (exact: both sides coincide)
    exact = True
    f_lo = f_hi = f
    pacc_lo = pacc_hi = pacc

    @property
    def threshold_bounds(self):
        return self.threshold_eps, self.threshold_eps


# ---------------------------------------------------------------------------
# Histograms (cached) -> curves
# ---------------------------------------------------------------------------
def _histograms(key, G0, H1, p):
    """(dim_A, A, cosets) for A = rowspan(G0) and its cosets A + h_i, cached in HIST_DIR."""
    path = HIST_DIR / f"{key}.json"
    if path.exists():
        h = json.loads(path.read_text())
        return h["dim_A"], h["A"], h["cosets"]
    G0 = np.asarray(G0, dtype=np.int64) % p
    H1 = np.asarray(H1, dtype=np.int64) % p
    dA, k = G0.shape[0], H1.shape[0]
    if p ** dA * max(1, k) <= LOCAL_NUMPY_WORDS:
        from .weightdist import coset_weight_enumerators

        hs = coset_weight_enumerators(G0, [None] + list(H1), p, max_words=int(LOCAL_NUMPY_WORDS))
        A, cosets = hs[0], hs[1:]
    else:
        import tempfile

        from dist_weightdist.coset_run import run_local

        t = max(1, min(dA - 1, int(math.ceil(math.log(64 * LOCAL_C_THREADS, p)))))
        with tempfile.TemporaryDirectory() as tmp:
            m = run_local(G0, H1, p, tmp, threads=LOCAL_C_THREADS, t=t)
        A, cosets = m["A"], m["cosets"]
    HIST_DIR.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"artifact_id": key, "p": p, "dim_A": dA,
                                "outputs": list(range(k)), "A": A, "cosets": cosets}))
    return dA, A, cosets


def _lambda_ok(H1, G0, p) -> bool:
    from .triorthogonal import is_cube_mod_p, logical_lambdas

    return all(is_cube_mod_p(l, p) for l in logical_lambdas(H1, len(H1), p))


def curve(name, G0, H1, p, key=None, note="") -> CodeCurve | None:
    """Curve for a code given (G0, H1); None if its logical gate is in another Clifford class."""
    H1 = np.atleast_2d(np.asarray(H1, dtype=np.int64) % p)
    if not _lambda_ok(H1, G0, p):
        return None
    dA, A, cosets = _histograms(key or name, G0, H1, p)
    prof = profile_from_histograms(p, dA, A, cosets)
    return CodeCurve(name=name, p=p, n=len(A) - 1, k=H1.shape[0], profile=prof, note=note)


def curve_from_record(record) -> CodeCurve | None:
    """Curve for a catalog record; uses a stored (possibly output-sampled) histogram if present."""
    path = HIST_DIR / f"{record.artifact_id}.json"
    note = record.family
    if path.exists():
        h = json.loads(path.read_text())
        prof = profile_from_histograms(record.p, h["dim_A"], h["A"], h["cosets"])
        if len(h["cosets"]) < record.k:
            note += f" (curve averaged over {len(h['cosets'])} sampled outputs)"
        c = CodeCurve(name=record.label, p=record.p, n=record.n, k=record.k, profile=prof, note=note)
        c.k = record.k
        return c
    from .performance import record_matrices

    G0, H1 = record_matrices(record)
    return curve(record.label, G0, H1, record.p, key=record.artifact_id, note=note)


# ---------------------------------------------------------------------------
# Library: catalog codes + the arXiv:2610.08012 cyclic / RS / RM codes
# ---------------------------------------------------------------------------
# (label, construction, length of the cyclic space, spectral support)  -- paper App. C
CYCLIC = {
    3: [("[[19,1,3]]_3 cyc-B", "B", 20, {2, 5, 6, 14, 15, 18}),
        ("[[40,1,4]]_3 cyc-A", "A", 40, {1, 2, 3, 5, 6, 9, 14, 15, 18, 27}),
        ("[[80,1,7]]_3 cyc-A", "A", 80, {1, 2, 3, 4, 5, 6, 9, 12, 15, 18, 27, 28, 36, 45, 54, 55}),
        ("[[79,1,7]]_3 cyc-B", "B", 80, {1, 3, 9, 14, 22, 27, 34, 38, 40, 42, 46, 53, 58, 66, 71, 77, 79}),
        ("[[121,1,8]]_3 cyc-A", "A", 121, {1, 2, 3, 6, 7, 9, 16, 18, 21, 23, 27, 41, 48, 54, 63, 68,
                                            69, 81, 83, 86})],
    5: [("[[11,1,3]]_5 cyc-B", "B", 12, {2, 3, 9, 10}),
        ("[[16,1,4]]_5 cyc-A", "A", 16, {1, 4, 5, 9, 13}),
        ("[[48,1,5]]_5 cyc-A", "A", 48, {1, 2, 3, 5, 10, 12, 15, 25, 27, 29, 39}),
        ("[[70,1,5]]_5 cyc-B", "B", 71, {1, 5, 14, 17, 25, 46, 54, 57, 66, 70}),
        ("[[90,1,6]]_5 cyc-B", "B", 91, {1, 5, 25, 28, 31, 34, 42, 47, 49, 51, 53, 63, 64, 73, 79, 83})],
    7: [("[[7,1,3]]_7 cyc-B", "B", 8, {1, 4, 7}),
        ("[[12,1,3]]_7 cyc-A", "A", 12, {1, 2, 7}),
        ("[[28,1,5]]_7 cyc-B", "B", 29, {1, 7, 16, 20, 23, 24, 25}),
        ("[[29,1,6]]_7 cyc-A", "A", 29, {1, 7, 16, 20, 23, 24, 25}),
        ("[[49,1,7]]_7 cyc-B", "B", 50, {1, 3, 7, 12, 16, 21, 25, 29, 34, 38, 43, 47, 49})],
}


def _rm_k1(p, m):
    from .performance import canonical_logical_rows
    from .reedmuller import r_max, rm_generator
    from .triorthogonal import build_triorthogonal_code

    r = r_max(m, p)
    G = rm_generator(r, m, p)
    G0 = build_triorthogonal_code(p, m, r, [1], G=G)["G0"]
    return G0, canonical_logical_rows(G, [0], p)


def library(p, include_catalog=True, verbose=False, dedupe=True) -> list[CodeCurve]:
    """All codes with an exact per-output curve for prime p (built once, histograms cached)."""
    from .cyclic import construction_A, construction_B

    out = []

    def add(c):
        if c is not None:
            out.append(c)
            if verbose:
                print(f"  {c.name:28s} n={c.n:4d} k={c.k:3d} {c.c_eps:.3g}e^{c.d} eps*={c.threshold_eps}",
                      flush=True)

    from .reedmuller import monomial_exponents, r_max

    for m in (1, 2) if p > 3 else (2, 3):   # shortest single-puncture RS / RM codes
        if p ** (len(monomial_exponents(m, r_max(m, p), p)) - 1) > LOCAL_NUMPY_WORDS:
            continue
        n = p ** m - 1
        tag = "RS" if m == 1 else "RM"
        add(curve(f"[[{n},1,?]]_{p} {tag}", *_rm_k1(p, m), p, key=f"rm-k1-p{p}-m{m}"))
    if p >= 7 and p not in CYCLIC:          # [[p,1,2*floor(p/6)+1]] family (paper Prop 15.2)
        nu, r = (p + 1) // 2, p // 6
        add(curve(f"[[{p},1,{2 * r + 1}]]_{p} cyc-B", *construction_B(p, p + 1, set(range(nu - r, nu + r + 1))),
                  p, key=f"cyclic-p{p}-B-n{p + 1}"))
    for label, cons, n, S in CYCLIC.get(p, []):
        G0, H1 = (construction_A if cons == "A" else construction_B)(p, n, S)
        add(curve(label, G0, H1, p, key=f"cyclic-p{p}-{cons}-n{n}"))
    if include_catalog:
        from .performance import assign_tier
        from .results import load_result_catalog

        seen = set()
        for r in load_result_catalog():
            if r.p != p or r.status != "confirmed" or r.m is None:
                continue
            sig = (r.n, r.k, r.A_d, r.puncture_columns)
            if sig in seen:
                continue
            if not (HIST_DIR / f"{r.artifact_id}.json").exists() and assign_tier(r)["tier"] != 1:
                continue
            seen.add(sig)
            add(curve_from_record(r))
    if not dedupe:
        return out
    # keep only the best curve per (n, k, d): near-duplicate codes only slow the search
    best = {}
    for c in out:
        key = (c.n, c.k, c.d)
        if key not in best or c.c_eps < best[key].c_eps:
            best[key] = c
    return list(best.values())


# ---------------------------------------------------------------------------
# Rigorous two-sided bounds for codes whose exact curve is out of reach
# ---------------------------------------------------------------------------
def _first_crossing(fun, top, grid=300, tol=1e-10):
    """Smallest eps in (0, top) with fun(eps) >= eps (log grid bracket + bisection), or None."""
    pts = list(np.logspace(-8, math.log10(top * (1 - 1e-6)), grid))
    prev = pts[0]
    if fun(prev) >= prev:
        return prev
    for e in pts[1:]:
        if fun(e) >= e:
            lo, hi = prev, e
            while hi - lo > tol * hi:
                mid = 0.5 * (lo + hi)
                lo, hi = (lo, mid) if fun(mid) >= mid else (mid, hi)
            return hi
        prev = e
    return None


class BoundedCurve:
    """Rigorous lower/upper bounds on one code's per-output map, from a subcode A' <= A.

    With d the exact distance (min weight of A^perp, so A^perp has no words of weight 1..d-1),
    A_d the number of weight-d logical vectors (None if unknown) and q = eps/(p-1):

      pacc_hi = P_acc(A')                       (A^perp <= A'^perp)
      pacc_lo = (1-eps)^n + A_d q^d (1-eps)^(n-d)
      f_lo    = (A_d/k) q^d (1-eps)^(n-d) / pacc_hi        (each weight-d logical hits >= 1 output;
                                                             A_d -> p-1 if unknown)
      f_hi    = min( eps_out(A')  [Lemma 17],
                     (A_d q^d (1-eps)^(n-d) + P(j in A'^perp, wt >= d+1)) / pacc_lo )
    The second f_hi term (needs A_d) is tight at small eps, where the subcode curve -- whose
    distance is smaller than d -- is loose.
    """

    exact = False

    def __init__(self, name, p, n, k, d, A_d, sub_profile: DistillationProfile, note=""):
        self.name, self.p, self.n, self.k, self.d, self.A_d = name, p, n, k, d, A_d
        self.sub = sub_profile
        self.note = note
        self._Blow = [sub_profile.accept.dual_count(w) for w in range(d + 1)]
        self.c_eps = (A_d if A_d else p - 1) / k / (p - 1) ** d   # lower-bound coefficient
        self.c_eps_hi = (A_d / (p - 1) ** d) if A_d else None
        top = (p - 1) / p
        self._thr = (_first_crossing(lambda e: float(self.f_hi(e)), top),
                     _first_crossing(lambda e: float(self.f_lo(e)), top))

    @property
    def threshold_bounds(self):
        return self._thr

    @property
    def threshold_eps(self):
        return self._thr[0]

    def _lead(self, eps, count):
        e = mpmath.mpf(eps)
        return count * (e / (self.p - 1)) ** self.d * (1 - e) ** (self.n - self.d)

    def pacc_hi(self, eps) -> float:
        with mpmath.workdps(self.sub._dps(float(eps))):
            return float(min(mpmath.mpf(1), self.sub.p_accept(eps)))

    def pacc_lo(self, eps) -> float:
        e = mpmath.mpf(eps)
        v = (1 - e) ** self.n + (self._lead(eps, self.A_d) if self.A_d else 0)
        return float(v)

    def f_lo(self, eps):
        with mpmath.workdps(self.sub._dps(float(eps))):
            return self._lead(eps, mpmath.mpf(self.A_d or self.p - 1) / self.k) / self.sub.p_accept(eps)

    def f_hi(self, eps):
        e = mpmath.mpf(eps)
        with mpmath.workdps(self.sub._dps(float(eps))):
            best = self.sub.eps_out(eps, "avg")
            if self.A_d:
                q = e / (self.p - 1)
                low = sum(self._Blow[w] * q ** w * (1 - e) ** (self.n - w) for w in range(self.d + 1))
                tail = max(mpmath.mpf(0), self.sub.p_accept(eps) - low)
                trunc = (self._lead(eps, self.A_d) + tail) / self.pacc_lo(eps)
                best = min(best, trunc)
            return best

    @property
    def gamma(self) -> float:
        return math.log(self.n / self.k) / math.log(self.d)


def bounded_curve_from_record(record, max_words=2e7, pool=12, seed=0) -> BoundedCurve:
    """Bounds for a catalog record: greedy Lemma-17 subcode (cached in HIST_DIR) + d, A_d."""
    from .performance import record_matrices, threshold_lower_bound_subcode

    assert record.d is not None, "needs the exact distance"
    path = HIST_DIR / f"bound-{record.artifact_id}.json"
    G0, H1 = record_matrices(record)
    if path.exists():
        rows = json.loads(path.read_text())["subcode"]
    else:
        b = threshold_lower_bound_subcode(G0, H1, record.p, kind="avg", max_words=max_words,
                                          pool=pool, seed=seed)
        rows = b["subcode"]
        HIST_DIR.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"artifact_id": record.artifact_id, "subcode": rows,
                                    "dim_A": int(G0.shape[0]), "max_words": max_words}))
    from .performance import profile_exact

    sub = profile_exact(np.array(rows), H1, record.p, kinds=("avg",), max_words=int(1e9))
    return BoundedCurve(record.label, record.p, record.n, record.k, record.d, record.A_d, sub,
                        note=f"bounds from subcode dim {len(rows)} of {G0.shape[0]}")


# ---------------------------------------------------------------------------
# Protocol search
# ---------------------------------------------------------------------------
@dataclass
class Protocol:
    codes: list
    cost: float
    eps: list   # eps after each round (eps[0] = input)

    def describe(self) -> str:
        steps = " -> ".join(c.name.split(" ")[0] for c in self.codes)
        return f"cost {self.cost:.4g}  [{steps}]  eps: " + ", ".join(mpmath.nstr(e, 3) for e in self.eps)


def best_protocol(curves, eps0, target, max_depth=6, only=None, mode="exact",
                  require=None) -> Protocol | None:
    """Cheapest sequence of codes taking per-output error eps0 to <= target.

    ``only`` restricts the search to codes whose name contains one of the given substrings.
    ``mode`` (for libraries containing ``BoundedCurve``s):
      * "upper": pessimistic curves (f_hi, pacc_lo) -> an upper bound on the optimal cost
        (the protocol found really reaches the target at no more than this cost);
      * "lower": optimistic curves (f_lo, pacc_hi) -> a lower bound on the optimal cost.
    Both assume each true eps_out is increasing in eps (checked on every exact curve);
    P_acc decreasing in eps is a theorem.
    ``require``: a curve that must appear in the sequence (cheapest protocol USING that code).
    """
    assert mode in ("exact", "upper", "lower")
    pool = [c for c in curves if only is None or any(s in c.name for s in only)]
    if mode == "exact":
        assert all(getattr(c, "exact", True) for c in pool), "bounded curves need mode upper/lower"
    fsel = {"exact": "f", "upper": "f_hi", "lower": "f_lo"}[mode]
    psel = {"exact": "pacc", "upper": "pacc_lo", "lower": "pacc_hi"}[mode]
    best = [None]

    def dfs(eps, cost, seq, traj):
        if eps <= target and (require is None or any(c is require for c in seq)):
            if best[0] is None or cost < best[0].cost:
                best[0] = Protocol(list(seq), cost, list(traj))
            return
        if len(seq) == max_depth:
            return
        for c in pool:
            pa = getattr(c, psel)(eps)
            if pa <= 0:
                continue
            new_cost = cost * c.n / (c.k * pa)
            if best[0] is not None and new_cost >= best[0].cost:
                continue
            e2 = getattr(c, fsel)(eps)
            if e2 >= eps:
                continue
            seq.append(c)
            traj.append(e2)
            dfs(e2, new_cost, seq, traj)
            seq.pop()
            traj.pop()

    dfs(mpmath.mpf(eps0), 1.0, [], [mpmath.mpf(eps0)])
    return best[0]


def protocol_bounds(curves, eps0, target, max_depth=6, only=None, require=None) -> dict:
    """Rigorous [lower, upper] bracket on the optimal protocol cost (see ``best_protocol``)."""
    lo = best_protocol(curves, eps0, target, max_depth, only, mode="lower", require=require)
    hi = best_protocol(curves, eps0, target, max_depth, only, mode="upper", require=require)
    return {"lower": lo, "upper": hi}


# ---------------------------------------------------------------------------
# One-code assessment: exact where possible, rigorous bounds otherwise
# ---------------------------------------------------------------------------
def code_curve(record, exact_words=3e8):
    """Exact curve if histograms are stored or the code is tier 1, else a BoundedCurve."""
    from .performance import assign_tier

    if (HIST_DIR / f"{record.artifact_id}.json").exists() or (
            assign_tier(record)["tier"] == 1
            and record.p ** assign_tier(record)["dim_A"] * max(1, record.k) <= exact_words):
        return curve_from_record(record)
    return bounded_curve_from_record(record)


def assess(record, eps0s=(1e-2, 1e-3), targets=(1e-10, 1e-15, 1e-30), max_depth=6) -> dict:
    """Leading order, threshold, P_acc and protocol cost of one catalog code.

    Cost = cheapest protocol that USES the code (other rounds drawn from ``library(p)``),
    compared with the cheapest protocol without it.  Exact when the code's curve is
    computable; otherwise every figure is a rigorous bound (see ``BoundedCurve``).
    """
    me = code_curve(record)
    def same(c):  # the code itself (an exact twin), not a different code with equal n, k, d
        return (c.n, c.k, c.d) == (me.n, me.k, me.d) and (
            not me.exact or abs(c.c_eps - me.c_eps) <= 1e-9 * me.c_eps)

    lib = [c for c in library(record.p, dedupe=False) if not same(c)]
    out = {"artifact_id": record.artifact_id, "label": record.label, "p": record.p,
           "exact": me.exact, "d": me.d, "note": me.note}
    if me.exact:
        out.update(c_eps=me.c_eps, threshold_eps=me.threshold_eps, P_acc_1e_3=me.pacc(1e-3))
    else:
        out.update(c_eps_bounds=(me.c_eps, me.c_eps_hi), threshold_eps_bounds=me.threshold_bounds,
                   P_acc_1e_3_bounds=(me.pacc_lo(1e-3), me.pacc_hi(1e-3)))
    rows = []
    for eps0 in eps0s:
        for tgt in targets:
            base = best_protocol(lib, eps0, tgt, max_depth)
            row = {"eps0": eps0, "target": tgt,
                   "without": base.cost if base else None,
                   "without_protocol": base.describe() if base else None}
            if me.exact:
                w = best_protocol(lib + [me], eps0, tgt, max_depth, require=me)
                row.update(with_cost=(w.cost, w.cost) if w else None,
                           with_protocol=w.describe() if w else None)
            else:
                b = protocol_bounds(lib + [me], eps0, tgt, max_depth, require=me)
                lo, hi = b["lower"], b["upper"]
                row.update(with_cost=(lo.cost if lo else None, hi.cost if hi else None),
                           with_protocol=hi.describe() if hi else None)
            rows.append(row)
    out["protocols"] = rows
    return out


def format_assessment(a) -> str:
    def rng(x, y, fmt="{:.4g}"):
        if x is None and y is None:
            return "n/a"
        if y is None:
            return ">= " + fmt.format(x)
        if x is not None and abs(y - x) <= 1e-6 * abs(y):
            return fmt.format(y)
        return f"[{fmt.format(x)}, {fmt.format(y)}]"

    lines = [f"{a['label']}  ({'exact' if a['exact'] else 'rigorous bounds: ' + a['note']})"]
    if a["exact"]:
        lines.append(f"  per-output eps_out ~ {a['c_eps']:.4g} eps^{a['d']}   threshold eps* = "
                     f"{a['threshold_eps']:.5g}   P_acc(1e-3) = {a['P_acc_1e_3']:.4g}")
    else:
        c_lo, c_hi = a["c_eps_bounds"]
        lines.append(f"  per-output eps_out ~ c eps^{a['d']}, c in {rng(c_lo, c_hi)}   threshold eps* "
                     f"{rng(*a['threshold_eps_bounds'], fmt='{:.4g}')}   P_acc(1e-3) in "
                     f"{rng(*a['P_acc_1e_3_bounds'])}")
    lines.append("  cost = raw magic states per output (cheapest protocol using this code vs without it)")
    for r in a["protocols"]:
        w = r["with_cost"]
        wtxt = rng(*w) if w else "unreachable"
        gain = ""
        if w and r["without"] and w[1] is not None:
            lo_c, hi_c, wo = w[0], w[1], r["without"]
            if wo > 1.005 * hi_c:
                tight = lo_c is not None and abs(hi_c - lo_c) <= 1e-6 * hi_c
                gain = f"   saving {'' if tight else '>= '}{wo / hi_c:.2f}x"
            elif lo_c is not None and wo > lo_c:
                gain = f"   undetermined (saving up to {wo / lo_c:.2f}x, or a loss)"
        lines.append(f"  {r['eps0']:g} -> {r['target']:g}:  with {wtxt}   without {r['without']:.4g}{gain}")
        if r["with_protocol"]:
            lines.append(f"      {r['with_protocol']}")
    return "\n".join(lines)
