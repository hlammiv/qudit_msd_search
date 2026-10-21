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


def library(p, include_catalog=True, verbose=False) -> list[CodeCurve]:
    """All codes with an exact per-output curve for prime p (built once, histograms cached)."""
    from .cyclic import construction_A, construction_B

    out = []

    def add(c):
        if c is not None:
            out.append(c)
            if verbose:
                print(f"  {c.name:28s} n={c.n:4d} k={c.k:3d} {c.c_eps:.3g}e^{c.d} eps*={c.threshold_eps}",
                      flush=True)

    for m in (1, 2) if p > 3 else (2, 3):   # shortest single-puncture RS / RM codes
        n = p ** m - 1
        tag = "RS" if m == 1 else "RM"
        add(curve(f"[[{n},1,?]]_{p} {tag}", *_rm_k1(p, m), p, key=f"rm-k1-p{p}-m{m}"))
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
    # keep only the best curve per (n, k, d): near-duplicate codes only slow the search
    best = {}
    for c in out:
        key = (c.n, c.k, c.d)
        if key not in best or c.c_eps < best[key].c_eps:
            best[key] = c
    return list(best.values())


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


def best_protocol(curves, eps0, target, max_depth=6, only=None) -> Protocol | None:
    """Cheapest sequence of codes taking per-output error eps0 to <= target.

    ``only`` restricts the search to codes whose name contains one of the given substrings.
    """
    pool = [c for c in curves if only is None or any(s in c.name for s in only)]
    best = [None]

    def dfs(eps, cost, seq, traj):
        if eps <= target:
            if best[0] is None or cost < best[0].cost:
                best[0] = Protocol(list(seq), cost, list(traj))
            return
        if len(seq) == max_depth:
            return
        for c in pool:
            pa = c.pacc(eps)
            if pa <= 0:
                continue
            new_cost = cost * c.n / (c.k * pa)
            if best[0] is not None and new_cost >= best[0].cost:
                continue
            e2 = c.f(eps)
            if e2 >= eps:
                continue
            seq.append(c)
            traj.append(e2)
            dfs(e2, new_cost, seq, traj)
            seq.pop()
            traj.pop()

    dfs(mpmath.mpf(eps0), 1.0, [], [mpmath.mpf(eps0)])
    return best[0]
