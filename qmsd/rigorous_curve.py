"""Rigorous per-output error bounds for a code whose exact curve cannot be enumerated.

For the large weighted-support codes (dim A ~ 56-66) the exact engine is out of reach and the
cost claims in PERFORMANCE.md use the leading-order surrogate c_eps eps^d.  This module turns a
*complete* weight-d census (A_d words, per-output average L_d) plus the certified distance into
rigorous two-sided bounds, so that a protocol found by ``best_protocol(..., mode="upper")``
provably reaches its target.

Noise model as in ``performance``: each input carries a nontrivial Z^j with probability eps,
uniformly over j != 0, so a specific error vector of weight w has probability
(eps/(p-1))^w (1-eps)^(n-w).  Accepted errors are the vectors of ker G0; output i fails on those
with h_i . e != 0.  With N_w the number of weight-w vectors of ker G0:

  * N_d = A_d exactly (census), and the per-output failing count is L_d on average;
  * for w > d, Singleton on the subcode supported inside any w-set (distance >= d, length w,
    hence dimension <= w-d+1) gives N_w <= C(n,w) (p^(w-d+1) - 1) < C(n,w) p^(w-d+1);
  * P_acc >= (1-eps)^n (the zero error) and P_acc <= (1-eps)^n + sum_{w>=d} N_w (...).

So  f_lo(eps) = L_d (eps/(p-1))^d (1-eps)^(n-d) / pacc_hi(eps)   (every weight-d failure counted,
P_acc at most pacc_hi) and f_hi(eps) = [L_d-term + sum_{w>d} C(n,w) p^(w-d+1) (eps/(p-1))^w
(1-eps)^(n-w)] / (1-eps)^n.  The tail is sum_{w>d} C(n,w) eta^w (1-eps)^(n-w) p^(1-d) with
eta = p eps/(p-1); it is summed exactly (mpmath) up to w = n.

The upper bound is loose by the Singleton factor (the random-code expectation of N_{d+1} is
C(n,d+1)(p-1)^(d+1)/p^(dim A), ~1e-50 of the bound), so it is 20+ orders above the surrogate at
eps ~ 1e-6 and forces the optimizer to run the big code at smaller input error.  That is the
price of a proof; the surrogate remains the realistic estimate.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache

import mpmath

DPS = 60


@dataclass
class RigorousBoundCurve:
    """Bounded curve (``exact = False``) for ``protocol.best_protocol(mode="upper"|"lower")``."""

    name: str
    p: int
    n: int
    k: int
    d: int
    L_d: float          # per-output average number of weight-d failing vectors (complete census)
    A_d: int            # total number of weight-d vectors of ker G0 (complete census)
    note: str = "rigorous bounds from a complete weight-d census + Singleton tail"
    # proved upper bounds for heavier weights, {w: (N_w bound, per-output failing bound)} for w = d+1, d+2, ...
    # (consecutive); the Singleton tail then starts after the last of them.  They enter f_hi and pacc_hi.
    heavier: dict = field(default_factory=dict)
    exact: bool = field(default=False, init=False)

    def __post_init__(self):
        ws = sorted(self.heavier)
        assert ws == list(range(self.d + 1, self.d + 1 + len(ws))), "heavier weights must be d+1, d+2, ... consecutive"
        for w, (N, L) in self.heavier.items():
            assert 0 <= L <= N, "per-output bound cannot exceed the word bound"

    @property
    def c_eps(self) -> float:
        return self.L_d / (self.p - 1) ** self.d

    @property
    def gamma(self) -> float:
        return math.log(self.n / self.k) / math.log(self.d)

    # ---- pieces
    def _vec(self, eps, w):
        e = mpmath.mpf(eps)
        return (e / (self.p - 1)) ** w * (1 - e) ** (self.n - w)

    def _known(self, eps, which):
        """sum over the proved heavier weights of bound * (eps/(p-1))^w (1-eps)^(n-w); which = 0 words, 1 per output."""
        return sum((mpmath.mpf(v[which]) * self._vec(eps, w) for w, v in self.heavier.items()), mpmath.mpf(0))

    def _tail(self, eps):
        """sum_{w>w0} C(n,w) p^(w-d+1) (eps/(p-1))^w (1-eps)^(n-w), unnormalised by P_acc, with w0 = d or the last
        weight in ``heavier``."""
        e = mpmath.mpf(eps); p, n, d = self.p, self.n, self.d
        w0 = d + len(self.heavier)
        eta = p * e / (p - 1); q = 1 - e
        acc = mpmath.mpf(0)
        if w0 >= n:
            return acc
        term = mpmath.binomial(n, w0 + 1) * eta ** (w0 + 1) * q ** (n - w0 - 1)
        for w in range(w0 + 1, n + 1):
            acc += term
            if w == n:
                break
            term = term * (n - w) / (w + 1) * eta / q
            if term < acc * mpmath.mpf(10) ** (-40) and w > n * float(eta):
                break
        return acc * mpmath.mpf(p) ** (1 - d)

    def _lead_vectors(self, eps):
        e = mpmath.mpf(eps)
        return (e / (self.p - 1)) ** self.d * (1 - e) ** (self.n - self.d)

    # ---- bounded interface
    def pacc_lo(self, eps) -> float:
        with mpmath.workdps(DPS):
            return float((1 - mpmath.mpf(eps)) ** self.n)

    def pacc_hi(self, eps) -> float:
        with mpmath.workdps(DPS):
            v = (1 - mpmath.mpf(eps)) ** self.n + self.A_d * self._lead_vectors(eps) + self._known(eps, 0) + self._tail(eps)
            return min(1.0, max(float(v), self.pacc_lo(eps)))

    def f_hi(self, eps):
        with mpmath.workdps(DPS):
            e = mpmath.mpf(eps)
            return (self.L_d * self._lead_vectors(eps) + self._known(eps, 1) + self._tail(eps)) / (1 - e) ** self.n

    def f_lo(self, eps):
        with mpmath.workdps(DPS):
            return self.L_d * self._lead_vectors(eps) / mpmath.mpf(self.pacc_hi(eps))

    def f_surrogate(self, eps):
        """The leading-order surrogate (objective.LeadingOrderCurve.f), for comparison."""
        e = mpmath.mpf(eps)
        return mpmath.mpf(self.c_eps) * (e / (1 - e)) ** self.d

    # conservative defaults for code that asks for a single curve
    f = f_hi
    pacc = pacc_lo

    @property
    def threshold_eps(self) -> float | None:
        """Largest eps with f_hi(eps) < eps (rigorous: below it one round provably improves)."""
        lo, hi = 1e-12, 0.5
        if self.f_hi(lo) >= lo:
            return None
        for _ in range(60):
            mid = math.sqrt(lo * hi)
            if self.f_hi(mid) < mid:
                lo = mid
            else:
                hi = mid
        return lo

    @property
    def threshold_bounds(self):
        return self.threshold_eps, None


def from_census(path_or_dict, name=None, heavier=None) -> RigorousBoundCurve:
    """Build from a census json (runs/weighted/*_census.json: code label, L_d, classes); ``heavier`` as in
    RigorousBoundCurve."""
    import json
    j = path_or_dict if isinstance(path_or_dict, dict) else json.load(open(path_or_dict))
    lab = j["code"]; n, k, d = [int(x) for x in lab.split("]]")[0][2:].split(",")]
    p = int(lab.split("_")[-1])
    A = j.get("A_d") or j.get("A_7") or j["classes"] * (p - 1)
    return RigorousBoundCurve(name or lab, p, n, k, d, float(j["L_d"]), int(A), heavier=dict(heavier or {}))


# ---------------------------------------------------------------------------------------------
# Weight-8 bound for the d = 7 sphere-charge codes at p = 7 (WEIGHTED_SUPPORTS.md, "Weight-8 words of RM(20,4)").
# A weight-8 kernel word is either a class-block word (<= p-1 per 8-set: the subcode on an 8-set has distance
# >= 7, so dimension <= 2, and a [8,2,7] MDS code has no weight-8 word), or a mixed word whose non-class part,
# reweighted by (g-1)(g-6), is a word of RM_7(20,4): of weight 8 (A_8 below, all conics) or a line word of weight
# 5-7 whose class completion of weight <= 3 is unique (two completions differ by a block word of weight <= 6).
RM20_A8 = 436_272_744_600                       # scratchpad/rm20_w8_count.py, runs/weighted/rm20_w8_count.json
RM20_A567 = 17_287_200 + 11_524_800 + 18_110_400


def sphere_weight8_bound(block_columns: int, p: int = 7) -> int:
    """Proved upper bound on the number of weight-8 kernel words of a d = 7 sphere-charge code (p = 7, m = 4,
    r = 7) whose class block has ``block_columns`` input columns."""
    assert p == 7
    return (p - 1) * math.comb(block_columns, 8) + RM20_A8 + RM20_A567


def _chain_ok(protocol, candidate, target) -> bool:
    """Independent re-check: feed eps through the chain with the candidate's f_hi and the exact
    curves of the others; the last eps must be <= target."""
    e = mpmath.mpf(protocol.eps[0])
    for c in protocol.codes:
        e = c.f_hi(e) if c is candidate else c.f(e)
    return e <= target


def rigorous_saving(candidate: RigorousBoundCurve, base, grid=None, max_depth=5) -> dict:
    """cost(base, exact) / cost(base + candidate) with the candidate's pessimistic (upper) and
    optimistic (lower) curves; the true optimal saving lies between ``saving_upper`` (proved)
    and ``saving_lower``."""
    from .objective import DEFAULT_GRID
    from .protocol import best_protocol
    rows = []
    for eps0, target in (grid or DEFAULT_GRID):
        b = best_protocol(base, eps0, target, max_depth=max_depth)
        up = best_protocol(base + [candidate], eps0, target, max_depth=max_depth, mode="upper")
        lo = best_protocol(base + [candidate], eps0, target, max_depth=max_depth, mode="lower")
        rows.append({"eps0": eps0, "target": target, "base": b and b.cost,
                     "with_upper": up and up.cost, "with_lower": lo and lo.cost,
                     "saving_proved": (b.cost / up.cost) if (b and up) else None,
                     "saving_optimistic": (b.cost / lo.cost) if (b and lo) else None,
                     "uses_proved": bool(up and any(c is candidate for c in up.codes)),
                     "protocol_proved": [c.name for c in up.codes] if up else None,
                     "eps_chain_proved": [float(x) for x in up.eps] if up else None,
                     "reaches_target_check": bool(up and _chain_ok(up, candidate, target))})
    return {"rows": rows}


# ---------------------------------------------------------------------------------------------
# Detect-and-correct rounds (2026-11-01).  Instead of discarding every block with a nonzero syndrome,
# decode to the lowest-weight error of weight <= t consistent with the syndrome (unique when 2t < d),
# undo it (a Z power is a Pauli, so the correction is free) and keep the block.
#   * acceptance >= P(wt(e) <= t)                                   (every such error is corrected)
#   * a kept block fails only if e = c + u with wt(c) <= t and u a nonzero kernel word with a nonzero
#     logical part, so wt(u) >= d.  Two rigorous bounds on P(fail and keep), the minimum is used:
#       (A) P(wt(e) >= d - t)
#       (B) sum_{w >= d} N_w sum_{s <= t} C(w, s) (eps/(p-1))^(w-s),  N_w <= C(n, w) p^(w-d+1) (Singleton)
#     ((B): e must agree with u on at least w - s of u's w support points.)
#   Per-output error is bounded by the block bound (conservative).  t = 0 reduces to pure detection with
#   the census-free Singleton bound.
@dataclass
class CorrectingBoundCurve:
    name: str
    p: int
    n: int
    k: int
    d: int
    t: int
    note: str = "detect-and-correct up to t errors; rigorous census-free bounds"
    exact: bool = field(default=False, init=False)

    def __post_init__(self):
        assert 2 * self.t < self.d, "unique decoding needs 2t < d"

    @property
    def gamma(self) -> float:
        return math.log(self.n / self.k) / math.log(self.d)

    def _binom_tail(self, eps, w0):
        e = mpmath.mpf(eps); n = self.n
        if w0 <= 0:
            return mpmath.mpf(1)
        term = mpmath.binomial(n, w0) * e ** w0 * (1 - e) ** (n - w0); acc = mpmath.mpf(0); w = w0
        while w <= n:
            acc += term
            if term < acc * mpmath.mpf(10) ** (-40) and w > n * float(e) + 5:
                break
            term = term * (n - w) / (w + 1) * e / (1 - e); w += 1
        return acc

    def _bound_B(self, eps):
        e = mpmath.mpf(eps); p, n, d, t = self.p, self.n, self.d, self.t; q = e / (p - 1)
        acc = mpmath.mpf(0); w = d
        while w <= n:
            inner = sum(mpmath.binomial(w, s) * q ** (w - s) for s in range(0, t + 1))
            term = mpmath.binomial(n, w) * mpmath.mpf(p) ** (w - d + 1) * inner
            acc += term
            if w > d + 3 and term < acc * mpmath.mpf(10) ** (-40):
                break
            w += 1
        return acc

    def pacc_lo(self, eps) -> float:
        with mpmath.workdps(DPS):
            return float(1 - self._binom_tail(eps, self.t + 1))

    def pacc_hi(self, eps) -> float:
        return 1.0

    def f_hi(self, eps):
        with mpmath.workdps(DPS):
            fail = min(self._binom_tail(eps, self.d - self.t), self._bound_B(eps))
            return fail / mpmath.mpf(self.pacc_lo(eps))

    def f_lo(self, eps):
        return mpmath.mpf(0)

    f = f_hi
    pacc = pacc_lo

    @property
    def threshold_eps(self):
        lo, hi = 1e-12, 0.5
        if self.f_hi(lo) >= lo:
            return None
        for _ in range(60):
            mid = math.sqrt(lo * hi)
            if self.f_hi(mid) < mid: lo = mid
            else: hi = mid
        return lo

    @property
    def threshold_bounds(self):
        return self.threshold_eps, None
