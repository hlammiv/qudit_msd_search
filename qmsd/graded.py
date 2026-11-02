"""Graded (detect-and-correct with recycling) distillation rounds, exact for small codes.

A round with code (G0, H1) at input error eps.  The exact joint distribution of (syndrome, logical
error) under the twirled channel is computed by a dynamic program over the n qudits on the state space
F_p^(r+k) (r = rank G0), so everything below is exact, not sampled or bounded:

  * the decoder maps each syndrome s to its lowest-weight explanation c(s) (coset leader, weight w(s));
  * a block is given *grade* w(s) if w(s) <= t, else discarded;
  * after the free Pauli correction, the block's logical error is (logical of e) - (logical of c(s)).

A graded round sends grades <= g forward, recycles grades g+1..t as inputs to the same round, and
discards the rest.  Recycling is modelled as a steady state: a fraction rho of each block's inputs are
recycled outputs (error y), the rest fresh (error x); the input error is the mixture (1-rho)x + rho y,
iterated to a fixed point.  This mixture step treats the recycled inputs as i.i.d. with the average
rate, which is exact for k = 1 codes in the limit of uniformly random placement and an approximation
otherwise; recycling is only offered for k = 1 codes.  Cost per forward output:

    fresh inputs per block / forward outputs per block = (n - k P_rec) / (k P_fwd),

exposed to ``protocol.best_protocol`` as an effective acceptance pacc_eff = n P_fwd / (n - k P_rec).
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

import numpy as np

from .structured_ad import _rank


def _independent_rows(G, p):
    G = np.asarray(G, dtype=np.int64) % p
    rows = []
    for row in G:
        if _rank(np.array(rows + [row]), p) > len(rows):
            rows.append(row)
    return np.array(rows, dtype=np.int64)


class GradedStats:
    """Exact grade statistics of one code; ``at(eps)`` returns per-grade acceptance and errors."""

    def __init__(self, name, G0, H1, p, tmax, max_states=2_000_000):
        self.name, self.p = name, p
        G0 = _independent_rows(G0, p); H1 = np.atleast_2d(np.asarray(H1, dtype=np.int64) % p)
        self.r, self.n = G0.shape; self.k = H1.shape[0]
        self.D = self.r + self.k
        self.size = p ** self.D
        assert self.size <= max_states, f"{name}: state space p^{self.D} = {self.size} too large"
        M = np.vstack([G0, H1]) % p                       # (r+k) x n
        self.cols = M.T.copy()                            # n x (r+k)
        self.base = p ** np.arange(self.D, dtype=np.int64)
        idx = np.arange(self.size, dtype=np.int64)
        self.digits = (idx[:, None] // self.base[None, :]) % p   # size x D
        # coset leaders up to weight tmax
        d_min = None
        leader = {0: (0, np.zeros(self.k, dtype=np.int64))}
        syn_base = self.base[: self.r]
        for w in range(1, tmax + 1):
            for S in itertools.combinations(range(self.n), w):
                for a in itertools.product(range(1, p), repeat=w):
                    v = (np.array(a) @ self.cols[list(S)]) % p
                    s = int(v[: self.r] @ syn_base)
                    if s == 0:
                        if v[self.r:].any() and d_min is None: d_min = w
                        continue
                    if s not in leader:
                        leader[s] = (w, v[self.r:].copy())
        self.leader = leader; self.tmax = tmax; self.d_detect = d_min; self._perms = {}
        nsyn = p ** self.r
        self._lw = -np.ones(nsyn, dtype=np.int64); self._ll = np.zeros((nsyn, self.k), dtype=np.int64)
        for s, (w, lc) in leader.items(): self._lw[s] = w; self._ll[s] = lc
        self._syn = self.digits[:, : self.r] @ self.base[: self.r]; self._logi = self.digits[:, self.r:]

    def _perm(self, j):
        if j not in self._perms:
            self._perms[j] = ((self.digits + self.cols[j][None, :]) % self.p) @ self.base
        return self._perms[j]

    def _distribution(self, eps):
        p = self.p; P = np.zeros(self.size); P[0] = 1.0
        q = eps / (p - 1)
        for j in range(self.n):
            perm = self._perm(j)            # state -> state + col_j (a bijection)
            new = (1 - eps) * P; cur = P
            for a in range(1, p):
                nxt = np.empty_like(cur); nxt[perm] = cur; cur = nxt     # shifted by a * col_j
                new += q * cur
            P = new
        return P

    def at(self, eps):
        """{grade: (probability, block logical-error prob, per-output logical-error prob)} at eps."""
        P = self._distribution(eps); p = self.p
        w = self._lw[self._syn]; ok = w >= 0
        resid = (self._logi - self._ll[self._syn]) % p
        bad_blk = (resid != 0).any(axis=1); bad_out = (resid != 0).sum(axis=1) / self.k
        out = {}
        for g in range(0, self.tmax + 1):
            m = ok & (w == g)
            if m.any():
                pr = P[m]; out[g] = (float(pr.sum()), float(pr[bad_blk[m]].sum()), float((pr * bad_out[m]).sum()))
        return out


@dataclass
class GradedCurve:
    """A graded round as a protocol curve (duck-typed like protocol.CodeCurve).

    t: correct up to weight t; g: grades <= g go forward; recycle: grades g+1..t re-enter as inputs."""
    stats: GradedStats
    t: int
    g: int
    recycle: bool
    grid_lo: float = 1e-9
    grid_hi: float = 0.3
    npts: int = 70
    note: str = "graded round (exact grade statistics; recycling as a steady-state mixture)"
    exact: bool = True
    _tab: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        s = self.stats
        self.p, self.n, self.k = s.p, s.n, s.k
        tag = f"t{self.t}g{self.g}" + ("r" if self.recycle else "")
        self.name = f"{s.name} [{tag}]"
        self.d = (s.d_detect or 0) - self.t if s.d_detect else 0
        xs = np.logspace(math.log10(self.grid_lo), math.log10(self.grid_hi), self.npts)
        fs, ps = [], []
        for x in xs:
            f, pa = self._solve(x); fs.append(max(f, 1e-300)); ps.append(max(pa, 1e-300))
        self._x = np.log(xs); self._f = np.log(fs); self._p = np.log(ps)
        self.threshold_eps = None
        for x, f in zip(xs[::-1], fs[::-1]):
            if f < x: self.threshold_eps = float(x); break

    def _round(self, eps):
        st = self.stats.at(eps)
        Pf = sum(v[0] for w, v in st.items() if w <= self.g)
        Ef = sum(v[2] for w, v in st.items() if w <= self.g)
        Pr = sum(v[0] for w, v in st.items() if self.g < w <= self.t) if self.recycle else 0.0
        Er = sum(v[2] for w, v in st.items() if self.g < w <= self.t) if self.recycle else 0.0
        return Pf, (Ef / Pf if Pf > 0 else 1.0), Pr, (Er / Pr if Pr > 0 else 1.0)

    def _solve(self, x):
        n, k = self.n, self.k
        if not self.recycle:
            Pf, f, _, _ = self._round(x); return f, Pf
        ein = x; rho = 0.0
        for _ in range(40):
            Pf, f, Pr, y = self._round(ein)
            if y >= x or Pr == 0:      # recycling would not help: fall back to discarding those grades
                Pf, f, _, _ = self._round(x); return f, Pf
            rho_new = min(k * Pr / n, 0.95)
            ein_new = (1 - rho_new) * x + rho_new * y
            if abs(ein_new - ein) < 1e-6 * x and abs(rho_new - rho) < 1e-9: break
            ein, rho = ein_new, rho_new
        pacc_eff = n * Pf / (n - k * Pr)
        return f, pacc_eff

    def _interp(self, arr, eps):
        le = math.log(eps)
        if le <= self._x[0]:   # extrapolate with the end slope (power law)
            sl = (arr[1] - arr[0]) / (self._x[1] - self._x[0]); return math.exp(arr[0] + sl * (le - self._x[0]))
        if le >= self._x[-1]: return math.exp(arr[-1])
        return math.exp(float(np.interp(le, self._x, arr)))

    def f(self, eps):
        import mpmath
        return mpmath.mpf(self._interp(self._f, float(eps)))

    def pacc(self, eps):
        return min(1.0, self._interp(self._p, float(eps)))

    f_lo = f_hi = f
    pacc_lo = pacc_hi = pacc

    @property
    def gamma(self):
        return math.log(self.n / self.k) / math.log(self.d) if self.d and self.d > 1 else float("inf")

    @property
    def threshold_bounds(self):
        return self.threshold_eps, self.threshold_eps
