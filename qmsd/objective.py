"""Protocol-cost search objective (Option C): how much a candidate code lowers distillation cost.

A search candidate is scored by its *marginal saving*: the cheapest multi-round protocol cost with
the standard library alone (RS / RM / cyclic / doubled / Prakash-Saha d=2), divided by the cost
when the candidate is added.  Exact curves are too expensive inside a search, so candidates use a
``LeadingOrderCurve`` built from (n, k, d) and the per-output weight-d count, which
``structured_ad`` gives cheaply:

    sum_i A_d^(i) = sum_w (w - d) * low_weight_histogram[w]

(each weight-w RM codeword with w - d points on the puncture set hits w - d outputs).

The surrogate keeps only the leading order, eps_out ~ c eps^d / (1 - eps)^d with P_acc ~
(1 - eps)^n, so it is optimistic at high input noise; ``validate`` compares it with exact curves.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import mpmath

from .protocol import STANDARD, best_protocol, library


@dataclass
class LeadingOrderCurve:
    """Duck-typed stand-in for ``protocol.CodeCurve`` from leading-order data only."""

    name: str
    p: int
    n: int
    k: int
    d: int
    L_d: float          # per-output average number of weight-d logical error vectors
    note: str = "leading-order surrogate"
    exact = True

    @property
    def c_eps(self) -> float:
        return self.L_d / (self.p - 1) ** self.d

    def f(self, eps):
        e = mpmath.mpf(eps)
        return mpmath.mpf(self.c_eps) * (e / (1 - e)) ** self.d

    def pacc(self, eps) -> float:
        return float((1 - mpmath.mpf(eps)) ** self.n)

    f_lo = f_hi = f
    pacc_lo = pacc_hi = pacc

    @property
    def threshold_eps(self) -> float | None:
        # c (e/(1-e))^d = e  ->  bisect on (0, (p-1)/p)
        g = lambda e: float(self.f(e)) - e
        lo, hi = 1e-12, (self.p - 1) / self.p - 1e-9
        if g(hi) < 0:
            return None
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            lo, hi = (mid, hi) if g(mid) < 0 else (lo, mid)
        return hi

    @property
    def gamma(self) -> float:
        return math.log(self.n / self.k) / math.log(self.d) if self.d > 1 else float("inf")


def per_output_Ld(low_weight_histogram: dict, d: int, k: int) -> float:
    """Per-output average weight-d count from ``structured_ad(...)['low_weight_histogram']``."""
    return sum((int(w) - d) * c for w, c in low_weight_histogram.items()) / k


def surrogate_from_structured_ad(p, m, r, puncture_columns, name=None) -> LeadingOrderCurve:
    from .structured_ad import structured_ad

    s = structured_ad(p, m, r, list(puncture_columns))
    k, d = len(puncture_columns), s["distance"]
    return LeadingOrderCurve(name or f"[[{p ** m - k},{k},{d}]]_{p}", p, p ** m - k, k, d,
                             per_output_Ld(s["low_weight_histogram"], d, k))


DEFAULT_GRID = ((1e-2, 1e-15), (1e-2, 1e-30), (1e-3, 1e-15), (1e-3, 1e-30))


def standard_library(p):
    return [c for c in library(p, include_catalog=False) if c.note == STANDARD]


def marginal_saving(candidate, base, grid=DEFAULT_GRID, max_depth=5) -> dict:
    """cost(base) / cost(base + candidate) at each (eps0, target); >= 1, == 1 if never used."""
    rows = []
    for eps0, target in grid:
        b = best_protocol(base, eps0, target, max_depth=max_depth)
        w = best_protocol(base + [candidate], eps0, target, max_depth=max_depth)
        ratio = (b.cost / w.cost) if (b and w) else (math.inf if w else None)
        rows.append({"eps0": eps0, "target": target, "base": b and b.cost, "with": w and w.cost,
                     "saving": ratio, "uses": bool(w and any(c is candidate for c in w.codes))})
    finite = [r["saving"] for r in rows if r["saving"] is not None]
    return {"best": max(finite) if finite else None, "rows": rows}


def validate(records, p, grid=DEFAULT_GRID, max_depth=5) -> list[dict]:
    """Surrogate vs exact marginal saving for catalog records that have exact curves."""
    from .performance import leading_order
    from .protocol import curve_from_record

    base = standard_library(p)
    out = []
    for r in records:
        exact = curve_from_record(r)
        lo = leading_order(exact.profile, "avg")
        sur = LeadingOrderCurve(r.label, r.p, r.n, r.k, lo.d, float(lo.L_d))
        out.append({"label": r.label, "exact": marginal_saving(exact, base, grid, max_depth)["best"],
                    "surrogate": marginal_saving(sur, base, grid, max_depth)["best"]})
    return out
