"""D = 3 -> 5 growth for Z_d colour-code cultivation (code deformation; Clifford proxy for T, zdsim).

The D=3 triangle sits in the corner of the D=5 triangle (identical stim coordinates).  Growth:
  1. prepare every new qudit in |0> (set Z) or |+> (set X), with the split chosen so that there are D=5
     logical representatives  Xbar = X^{xbar}, Zbar = Z^{zbar}  with  xbar|_old ~ X_L^(3), zbar|_old ~ Z_L^(3)
     (mod D=3 stabilisers), xbar = 0 on Z-set and zbar = 0 on X-set.  Then Xbar, Zbar act on the prepared
     state exactly as the D=3 logicals, and they commute with every D=5 stabiliser, so measuring the D=5
     stabilisers preserves the encoded (magic) state;
  2. measure all D=5 plaquettes; outcomes that are deterministic get detectors, random ones are
     gauge-fixed by a feed-forward Pauli (destabiliser) that flips only that plaquette and commutes with
     Xbar and Zbar (zdsim op FF);
  3. continue with D=5 stabiliser rounds and T-checks; ideal D=5 decoding at the end.
Logical preservation is verified operationally: every later detector (in particular the D=5 T-checks) must be
deterministic (zdsim.check_determinism).
"""

from __future__ import annotations

import itertools as it

import numpy as np

from .cult_circuit import ColorCult, pivot_form
from .cultivation import color_code_plaquettes, honeycomb_signs
from .ringlinalg import right_kernel_howell
from .zdsim import Circuit, FrameSim, check_determinism


def embedding(D_small: int = 3, D_big: int = 5):
    """indices of the small code's qudits inside the big one, matched by stim coordinates."""
    Hs, xs = color_code_plaquettes(D_small)
    Hb, xb = color_code_plaquettes(D_big)
    look = {tuple(np.round(p, 6)): i for i, p in enumerate(xb)}
    emb = [look[tuple(np.round(p, 6))] for p in xs]
    return Hs, Hb, xb, emb


def _solve_affine(A, b, d):
    """some u with A u = b (mod d), or None (kernel of [A | -b] with a unit last coordinate)."""
    Aug = np.concatenate([A % d, (-np.asarray(b)[:, None]) % d], axis=1)
    K = right_kernel_howell(Aug, d)
    for v in K:
        if np.gcd(int(v[-1]) % d, d) == 1:
            v = np.array(v, dtype=np.int64)
            return (v[:-1] * pow(int(v[-1]) % d, -1, d)) % d
    return None


def z_logical(H, MZ, d):
    """a Z-type logical c (H c = 0) with c . 1 = 1 (pairs with X_L = X^{all ones})."""
    n = H.shape[1]
    A = np.concatenate([H, np.ones((1, n), dtype=np.int64)])
    b = np.zeros(A.shape[0], dtype=np.int64); b[-1] = 1
    return _solve_affine(A, b, d)


def logical_reps(Hb, MZb, Hs, MZs, emb, zset, d):
    """(xbar, zbar) for a new-qudit basis split (zset = new qudits prepared in |0>), or None."""
    nb, ns = Hb.shape[1], Hs.shape[1]
    new = [i for i in range(nb) if i not in emb]
    xset = [i for i in new if i not in zset]
    # xbar: unknowns (xbar in Z_d^nb, w in Z_d^{rows Hs});  MZb xbar = 0, xbar|old - w Hs = 1, xbar[zset] = 0
    rs = Hs.shape[0]
    rows, rhs = [], []
    for r in MZb:
        rows.append(np.concatenate([r, np.zeros(rs, dtype=np.int64)])); rhs.append(0)
    for j, i in enumerate(emb):
        e = np.zeros(nb + rs, dtype=np.int64); e[i] = 1; e[nb:] = (-Hs[:, j]) % d
        rows.append(e); rhs.append(1)
    for i in zset:
        e = np.zeros(nb + rs, dtype=np.int64); e[i] = 1; rows.append(e); rhs.append(0)
    u = _solve_affine(np.array(rows), np.array(rhs), d)
    if u is None:
        return None
    xbar = u[:nb]
    zl = z_logical(Hs, MZs, d)
    rows, rhs = [], []
    for r in Hb:
        rows.append(np.concatenate([r, np.zeros(MZs.shape[0], dtype=np.int64)])); rhs.append(0)
    for j, i in enumerate(emb):
        e = np.zeros(nb + MZs.shape[0], dtype=np.int64); e[i] = 1; e[nb:] = (-MZs[:, j]) % d
        rows.append(e); rhs.append(int(zl[j]))
    for i in xset:
        e = np.zeros(nb + MZs.shape[0], dtype=np.int64); e[i] = 1; rows.append(e); rhs.append(0)
    u = _solve_affine(np.array(rows), np.array(rhs), d)
    if u is None:
        return None
    return xbar, u[:nb]


def destabilizers(Hrows, other_logical, d):
    """E[s]: vector with Hrows E[s] = delta_s and other_logical . E[s] = 0 (mod d)."""
    A = np.concatenate([Hrows, other_logical[None, :]])
    out = []
    for s in range(Hrows.shape[0]):
        b = np.zeros(A.shape[0], dtype=np.int64); b[s] = 1
        e = _solve_affine(A, b, d)
        if e is None:
            raise ValueError("no destabiliser")
        out.append(e)
    return np.array(out)


class GrowCult(ColorCult):
    """Inject + cultivate at D=3, grow to D=5, cultivate at D=5.

    rounds3/checks3: D=3 stage (stabiliser round + checks, repeated); rounds5/checks5 likewise after growth.
    zset: new qudits prepared in |0> (default: first feasible split found by `growth_splits`, preferring the
    most deterministic plaquettes)."""

    def __init__(self, d: int = 4, p: float = 1e-3, sigma: int = 1, inj: int | None = None,
                 rounds3: int = 1, checks3: int = 2, rounds5: int = 1, checks5: int = 2,
                 zset=None, noisy_injection: bool = True):
        Hs, Hb, xb, emb = embedding(3, 5)
        tau5 = honeycomb_signs(Hb, xb)
        self.d, self.p, self.sigma = d, p, sigma
        self.inj = sigma if inj is None else inj
        self.emb = emb
        self.Hs, self.Hb, self.tau5 = Hs, Hb, tau5
        tau3 = tau5[emb]
        assert not np.any(Hs @ tau3), "D=3 restriction of the D=5 bipartition is unbalanced"
        self.MZs, self.MZb = (Hs * tau3) % d, (Hb * tau5) % d
        if zset is None:
            zset = growth_splits(d, limit=1)[0]["zset"]
        self.zset = list(zset)
        reps = logical_reps(Hb, self.MZb, Hs, self.MZs, emb, self.zset, d)
        if reps is None:
            raise ValueError("infeasible split")
        self.xbar, self.zbar = reps
        c = self.c = Circuit(d)
        self.data5 = c.q(Hb.shape[1])
        xa5, za5, cat5 = c.q(len(Hb)), c.q(len(Hb)), c.q(Hb.shape[1])
        # ---- D=3 stage
        self._phase(Hs, tau3, [self.data5[i] for i in emb], xa5[:len(Hs)], za5[:len(Hs)], cat5[:len(emb)], 3)
        self.encode(noisy_injection)
        for _ in range(rounds3):
            self.stab_round()
            for _ in range(checks3):
                self.check()
        # ---- growth
        self._phase(Hb, tau5, self.data5, xa5, za5, cat5, 5)
        self.grow()
        for _ in range(rounds5):
            self.stab_round()
            for _ in range(checks5):
                self.check()
        self.stab_round()

    def _phase(self, H, tau, data, xa, za, cat, D):
        self.H, self.tau, self.data, self.xa, self.za, self.cat = H, tau, list(data), list(xa), list(za), list(cat)
        self.n, self.distance, self.MZ, self._tab = H.shape[1], D, (H * tau) % self.d, None

    def grow(self):
        d, H, MZ, D = self.d, self.H, self.MZ, self.data
        new = [i for i in range(self.n) if i not in self.emb]
        self.R([D[i] for i in new])
        xs = [D[i] for i in new if i not in self.zset]
        if xs:
            self.F(xs)
        # first D=5 round: measure, then classify deterministic vs random in a noiseless gauge run
        m_x, m_z = [], []
        for r in range(len(H)):
            a, b = self.xa[r], self.za[r]
            sup = list(np.nonzero(H[r])[0])
            self.R([a, b]); self.F([a])
            self.CSUM([(a, D[i]) for i in sup], 1)
            self.CSUM([(D[i], b) for i in sup], [int(self.tau[i]) % d for i in sup])
            self.FD([a])
            ma, mb = self.M([a, b])
            m_x.append(ma); m_z.append(mb)
        s = FrameSim(self.c, 128, np.random.default_rng(7), gauge=True).run(scale=0.0)
        rnd = np.any(s.meas, axis=0)
        EZ = destabilizers(H, self.xbar, d)          # Z-type fixes for X-plaquettes (commute with Xbar)
        EX = destabilizers(MZ, self.zbar, d)         # X-type fixes for Z-plaquettes (commute with Zbar)
        ff = []
        self.random_plaquettes = {"X": [], "Z": []}
        for r in range(len(H)):
            if rnd[m_x[r]]:
                self.random_plaquettes["X"].append(r)
                ff += [(m_x[r], D[q], "z", int(EZ[r, q])) for q in range(self.n) if EZ[r, q]]
            else:
                self.c.detector([(m_x[r], 1)], "grow-X")
            if rnd[m_z[r]]:
                self.random_plaquettes["Z"].append(r)
                ff += [(m_z[r], D[q], "x", int(-EX[r, q]) % d) for q in range(self.n) if EX[r, q]]
            else:
                self.c.detector([(m_z[r], 1)], "grow-Z")
        self.c.FF(ff)


def growth_splits(d: int = 4, limit: int | None = None):
    """Feasible |0>/|+> splits of the 12 new qudits, sorted by the number of random (gauge) plaquettes."""
    Hs, Hb, xb, emb = embedding(3, 5)
    tau5 = honeycomb_signs(Hb, xb)
    MZs, MZb = (Hs * tau5[emb]) % d, (Hb * tau5) % d
    new = [i for i in range(Hb.shape[1]) if i not in emb]
    feas = []
    for bits in it.product((0, 1), repeat=len(new)):
        zset = [q for q, b in zip(new, bits) if b]
        if logical_reps(Hb, MZb, Hs, MZs, emb, zset, d) is None:
            continue
        # count random plaquettes: X-plaquette deterministic iff in span(old X stabs, X on xset); Z likewise
        xset = [q for q in new if q not in zset]
        nr = 0
        for r in range(len(Hb)):
            for rows, row in ((_span_rows(Hs, emb, xset, Hb.shape[1], d), Hb[r]),
                              (_span_rows(MZs, emb, zset, Hb.shape[1], d), MZb[r])):
                if _solve_affine(rows.T, row, d) is None:
                    nr += 1
        feas.append({"zset": zset, "random_plaquettes": nr})
    feas.sort(key=lambda f: f["random_plaquettes"])
    return feas[:limit] if limit else feas


def _span_rows(Hs, emb, singles, nb, d):
    rows = []
    for r in Hs:
        v = np.zeros(nb, dtype=np.int64); v[emb] = r; rows.append(v)
    for q in singles:
        v = np.zeros(nb, dtype=np.int64); v[q] = 1; rows.append(v)
    return np.array(rows)
