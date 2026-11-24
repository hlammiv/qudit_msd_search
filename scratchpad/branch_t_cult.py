"""Branch-sum exact-T simulator for Z_d cultivation (scalable successor of exact_t_cult's dense data register).

Exact for Pauli noise.  As in exact_t_cult, a Pauli frame carries every Clifford part and the data frame is flushed at
segment boundaries; the data register is represented not as a state vector but as a short sum of Pauli branches on the
IDEAL (noiseless) data state psi(t):

    state = sum_i a_i P_i psi(t),   P_i = X^x Z^z (X left of Z, w = e^{2 pi i/d}, Z X = w X Z).

  frame flush F : P_i -> F P_i                       (phase w^{z_F . x_i})
  T layer T^t   : T X^x T^-1 = X^x sum_m c_m(x, t) Z^m  per qudit (branches multiply by d per X-hit qudit); psi -> T psi
  detector m    : P_i psi has the definite eigen-index shift v.x - u.z of O_m = X^u Z^v (psi is an eigenstate), so the
                  post-selection keeps the rows with shift + frame flip = 0 (no amplitude change)
  data Clifford : exact conjugation tables with phases (dense d x d / d^2 x d^2 matrices)
  random meas.  : (GrowCult's growth round, gauge-fixed by feed-forward E_m) -- with G_m = Z^{g_m} on the fresh qudits
                  (a stabiliser of psi flipping only m), the sum over outcome records is exact and gives:
                  row P_i -> P_i G^{-shift_i}, and rows with different s_{i,m} = <E_m, P_i G^{-shift_i}> stop
                  interfering (incoherent classes, `cls`).
  output        : rows grouped by final syndrome sector; sector rep R_s = X^e Z^e' (decoder corrections);
                  P_i psi = R_s w^{ph_i} X_L^a Z_L^b psi; with beta_ab the summed amplitudes:
                  w_s = beta^+ G beta (Gram of the logical state), F_s = |sum beta_ab e(a,b)|^2,
                  e(a, b) = <psi| X^{a xl} Z^{b zl} |psi> by exact backward conjugation through the encoder.
Validation (`selftest`): ColorCult D=3 (d = 2, 4) and D=5 (d = 2) shot by shot against exact_t_cult's dense hybrid;
GrowCult noiseless (w = F = 1) and in proxy mode (T -> S) shot by shot against the Clifford frame proxy.

    PYTHONPATH=.:scratchpad python scratchpad/branch_t_cult.py selftest
    PYTHONPATH=.:scratchpad python scratchpad/branch_t_cult.py run d p shots out.json seed [batch] [cfg]
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

from primepower_msd.cult_circuit import ColorCult
from primepower_msd.cult_growth import GrowCult, _solve_affine
from primepower_msd.cultivation import escaping_carrier
from primepower_msd.ringlinalg import howell_form, in_span, right_kernel_howell
from primepower_msd.zdsim import Circuit, FrameSim, add_tick_idle


# ----------------------------------------------------------------------------------------------- exact Pauli tables
class Tables:
    def __init__(self, d, q, proxy=False):
        self.d, self.q, self.proxy = d, np.asarray(q), proxy
        self.w = np.exp(2j * np.pi / d)
        self.Xm = np.roll(np.eye(d), 1, axis=0)
        self.Zm = np.diag(self.w ** np.arange(d))
        self.Fm = np.array([[self.w ** (j * k) for j in range(d)] for k in range(d)]) / np.sqrt(d)
        self._c = {}
        self._t = {}

    def P1(self, x, z):
        return np.linalg.matrix_power(self.Xm, x % self.d) @ np.linalg.matrix_power(self.Zm, z % self.d)

    def gate(self, key):
        d = self.d
        if key == "F":
            return self.Fm
        if key == "FD":
            return self.Fm.conj().T
        if key[0] == "S":
            x = np.arange(d)
            return np.diag(np.exp(2j * np.pi * ((key[1] * x * x) % (2 * d)) / (2 * d)))
        if key[0] == "CSUM":
            w = key[1]
            G = np.zeros((d * d, d * d))
            for a in range(d):
                for b in range(d):
                    G[a * d + (b + w * a) % d, a * d + b] = 1
            return G
        raise KeyError(key)

    def conj(self, key):
        """table for P -> G P G^dag: arrays (tx, tz, tc) over flat Pauli index (x z per qudit, 1 or 2 qudits)."""
        if key in self._c:
            return self._c[key]
        d = self.d
        G = self.gate(key)
        nq = 1 if G.shape[0] == d else 2
        N = d ** (2 * nq)
        tx = np.zeros((N, nq), np.int16); tz = np.zeros((N, nq), np.int16); tc = np.zeros(N, complex)
        for idx in range(N):
            if nq == 1:
                x, z = divmod(idx, d); P = self.P1(x, z)
            else:
                x1, z1, x2, z2 = (idx // d ** 3) % d, (idx // d ** 2) % d, (idx // d) % d, idx % d
                P = np.kron(self.P1(x1, z1), self.P1(x2, z2))
            M = G @ P @ G.conj().T
            col = M[:, 0]
            r = int(np.argmax(np.abs(col)))
            c0 = col[r]
            if nq == 1:
                xs = [r]
                z1_ = int(np.round(np.angle(M[(1 + r) % d, 1] / c0) / (2 * np.pi / d))) % d
                zs = [z1_]
                Q = self.P1(xs[0], zs[0])
            else:
                xs = [r // d, r % d]
                r10 = ((1 + xs[0]) % d) * d + xs[1]; r01 = xs[0] * d + (1 + xs[1]) % d
                zs = [int(np.round(np.angle(M[r10, d] / c0) / (2 * np.pi / d))) % d,
                      int(np.round(np.angle(M[r01, 1] / c0) / (2 * np.pi / d))) % d]
                Q = np.kron(self.P1(xs[0], zs[0]), self.P1(xs[1], zs[1]))
            c = c0 / Q[r, 0]
            assert np.allclose(M, c * Q, atol=1e-9), (key, idx)
            tx[idx] = xs; tz[idx] = zs; tc[idx] = c
        self._c[key] = (tx, tz, tc)
        return self._c[key]

    def tcoef(self, t):
        """C[x, m]: T^t X^x T^-t = X^x sum_m C[x, m] Z^m (exact T, or S^t in proxy mode)."""
        if t in self._t:
            return self._t[t]
        d = self.d
        y = np.arange(d)
        C = np.zeros((d, d), complex)
        for x in range(d):
            if self.proxy:
                ph = lambda v: np.exp(2j * np.pi * ((t * v * v) % (2 * d)) / (2 * d))
            else:
                ph = lambda v: np.exp(2j * np.pi * ((t * self.q[v % d]) % (4 * d)) / (4 * d))
            D = ph((y + x) % d) / ph(y)
            C[x] = np.array([np.sum(D * self.w ** (-m * y)) for m in range(d)]) / d
        C[np.abs(C) < 1e-13] = 0
        self._t[t] = C
        return C


# ----------------------------------------------------------------------------------------------- rows (branches)
class Rows:
    def __init__(self, n, sids, d):
        R = len(sids)
        self.d = d
        self.X = np.zeros((R, n), np.int16); self.Z = np.zeros((R, n), np.int16)
        self.a = np.ones(R, complex); self.sid = np.asarray(sids, np.int64); self.cls = np.zeros(R, np.int64)
        self.shift = {}

    def take(self, idx):
        self.X, self.Z, self.a, self.sid, self.cls = self.X[idx], self.Z[idx], self.a[idx], self.sid[idx], self.cls[idx]
        self.shift = {m: s[idx] for m, s in self.shift.items()}

    def __len__(self):
        return len(self.a)

    def left(self, xF, zF):
        """rows <- F P (xF, zF per row)."""
        d = self.d
        self.a = self.a * np.exp(2j * np.pi * ((np.einsum("ij,ij->i", zF.astype(np.int64), self.X.astype(np.int64))) % d) / d)
        self.X = ((self.X + xF) % d).astype(np.int16); self.Z = ((self.Z + zF) % d).astype(np.int16)

    def right(self, xG, zG):
        """rows <- P G."""
        d = self.d
        self.a = self.a * np.exp(2j * np.pi * ((np.einsum("ij,ij->i", self.Z.astype(np.int64), xG.astype(np.int64))) % d) / d)
        self.X = ((self.X + xG) % d).astype(np.int16); self.Z = ((self.Z + zG) % d).astype(np.int16)

    def conj1(self, j, tab):
        tx, tz, tc = tab
        idx = self.X[:, j].astype(np.int64) * self.d + self.Z[:, j]
        self.X[:, j] = tx[idx, 0]; self.Z[:, j] = tz[idx, 0]; self.a = self.a * tc[idx]

    def conj2(self, j, k, tab):
        tx, tz, tc = tab
        d = self.d
        idx = ((self.X[:, j].astype(np.int64) * d + self.Z[:, j]) * d + self.X[:, k]) * d + self.Z[:, k]
        self.X[:, j] = tx[idx, 0]; self.Z[:, j] = tz[idx, 0]; self.X[:, k] = tx[idx, 1]; self.Z[:, k] = tz[idx, 1]
        self.a = self.a * tc[idx]

    def texpand(self, j, C):
        hit = np.nonzero(self.X[:, j] != 0)[0]
        if not len(hit):
            return
        keep = np.nonzero(self.X[:, j] == 0)[0]
        xs = self.X[hit, j].astype(np.int64)
        parts = [keep]
        news = []
        d = self.d
        coef = C[xs]                                   # (H, d)
        ii, mm = np.nonzero(coef)
        src = hit[ii]
        old = len(self.a)
        self.take(np.concatenate([keep, src]))
        nk = len(keep)
        self.Z[nk:, j] = (self.Z[nk:, j] + mm) % d
        self.a[nk:] = self.a[nk:] * coef[ii, mm]

    def reduce_z(self, red):
        """Z part modulo a Z-stabiliser group of the ideal state (Howell rows, eigenvalue 1, no phase)."""
        if red is None or not len(self.a):
            return
        d = self.d
        Z = self.Z.astype(np.int64)
        for c, pv, row in red:
            q = Z[:, c] // pv
            if q.any():
                Z = (Z - q[:, None] * row[None, :]) % d
        self.Z = Z.astype(np.int16)

    def reduce_x(self, gens, sweeps=3):
        """lower the X weight of each row by X-stabilisers h of the ideal state (X^h psi = psi):
        X^{x+u} Z^z psi = w^{-u.z} X^x Z^z psi.  Greedy weight descent over the generators and their multiples."""
        if gens is None or not len(self.a):
            return
        d = self.d
        X = self.X.astype(np.int64); Z = self.Z.astype(np.int64)
        U = np.zeros_like(X)
        for _ in range(sweeps):
            changed = False
            for h in gens:
                for c in range(1, d):
                    cand = (X - c * h[None, :]) % d
                    better = (cand != 0).sum(1) < (X != 0).sum(1)
                    if better.any():
                        X[better] = cand[better]; U[better] = (U[better] + c * h[None, :]) % d; changed = True
            if not changed:
                break
        self.a = self.a * np.exp(2j * np.pi * ((-(U * Z).sum(1)) % d) / d)
        self.X = X.astype(np.int16)

    def merge(self, red=None):
        self.reduce_z(red)
        if len(self.a) < 2:
            return
        key = np.concatenate([self.sid[:, None].view(np.uint8).reshape(len(self.a), -1),
                              self.cls[:, None].view(np.uint8).reshape(len(self.a), -1),
                              self.X.view(np.uint8), self.Z.view(np.uint8)], axis=1)
        key = np.ascontiguousarray(key).view(np.dtype((np.void, key.shape[1]))).ravel()
        u, first, inv = np.unique(key, return_index=True, return_inverse=True)
        if len(u) == len(key):
            return
        ar = np.bincount(inv, weights=self.a.real, minlength=len(u)); ai = np.bincount(inv, weights=self.a.imag, minlength=len(u))
        self.take(first)
        self.a = ar + 1j * ai
        nz = np.abs(self.a) > 1e-13
        if not nz.all():
            self.take(np.nonzero(nz)[0])


# ----------------------------------------------------------------------------------------------- frame pass (with FF)
def frame_pass(c: Circuit, data, flush_ops, B, rng, scale=1.0, forced=None, skip_ff=False):
    """exact_t_cult.frame_pass + feed-forward (FF on the recorded flips of each frame; skip_ff for observable runs) and
    data resets anywhere (a flush op that is a data R: recorded before the reset)."""
    d = c.d
    x = np.zeros((B, c.nq), np.int64); z = np.zeros_like(x)
    xp = np.zeros_like(x); zp = np.zeros_like(x)
    meas = np.zeros((B, c.n_meas), np.int64); measp = np.zeros_like(meas)
    rec = {}
    data = list(data)
    for k, op in enumerate(c.ops):
        name = op[0]
        if k in flush_ops:
            rec[k] = (x[:, data].copy(), z[:, data].copy())
            x[:, data] = 0; z[:, data] = 0
        for X_, Z_, M_ in ((x, z, meas), (xp, zp, measp)):
            if name == "R":
                X_[:, op[1]] = 0; Z_[:, op[1]] = 0
            elif name == "F":
                q = op[1]; X_[:, q], Z_[:, q] = (-Z_[:, q]) % d, X_[:, q].copy()
            elif name == "FD":
                q = op[1]; X_[:, q], Z_[:, q] = Z_[:, q].copy(), (-X_[:, q]) % d
            elif name == "S":
                for qq, a in zip(op[1], op[2]):
                    Z_[:, qq] = (Z_[:, qq] + a * X_[:, qq]) % d
            elif name == "CSUM":
                for (cq, tq), w in zip(op[1], op[2]):
                    X_[:, tq] = (X_[:, tq] + w * X_[:, cq]) % d
                    Z_[:, cq] = (Z_[:, cq] - w * Z_[:, tq]) % d
            elif name == "M":
                M_[:, op[2]] = X_[:, op[1]]
            elif name == "FF" and not skip_ff:
                for i, qq, ax, cf in op[1]:
                    arr = X_ if ax == "x" else Z_
                    arr[:, qq] = (arr[:, qq] + cf * M_[:, i]) % d
        if name in ("DEP1", "XERR", "ZERR", "DEP2") and op[-1] * scale > 0:
            p = op[-1] * scale
            if name == "DEP2":
                pairs = op[1]
                hit = rng.random((B, len(pairs))) < p
                r = rng.integers(1, d ** 4, size=hit.shape)
                P = np.stack([(r // d ** j) % d for j in range(4)], axis=-1)
                a = [pp[0] for pp in pairs]; b = [pp[1] for pp in pairs]
                for X_, Z_ in ((x, z), (xp, zp)):
                    X_[:, a] = (X_[:, a] + hit * P[..., 0]) % d; Z_[:, a] = (Z_[:, a] + hit * P[..., 1]) % d
                    X_[:, b] = (X_[:, b] + hit * P[..., 2]) % d; Z_[:, b] = (Z_[:, b] + hit * P[..., 3]) % d
            else:
                q = op[1]
                hit = rng.random((B, len(q))) < p
                if name == "DEP1":
                    r = rng.integers(1, d ** 2, size=hit.shape)
                    P0, P1 = r % d, (r // d) % d
                else:
                    v = rng.integers(1, d, size=hit.shape)
                    P0, P1 = (v, 0 * v) if name == "XERR" else (0 * v, v)
                for X_, Z_ in ((x, z), (xp, zp)):
                    X_[:, q] = (X_[:, q] + hit * P0) % d; Z_[:, q] = (Z_[:, q] + hit * P1) % d
        if forced is not None and k in forced:
            for mask, qs, P in forced[k]:
                for X_, Z_ in ((x, z), (xp, zp)):
                    X_[:, qs] = (X_[:, qs] + mask[:, None] * P[..., 0]) % d
                    Z_[:, qs] = (Z_[:, qs] + mask[:, None] * P[..., 1]) % d
    rec["end"] = (x[:, data].copy(), z[:, data].copy())
    return dict(rec=rec, flips=meas % d, xp=xp[:, data], zp=zp[:, data], measp=measp % d)


# ----------------------------------------------------------------------------------------------- plan
class Plan:
    """c, data (all data qudits, final code order), tlayers {op: {qudit: t}}, carrier q; final: dict(H, MZ, xl, zl,
    tables=(tx, tz)); e_fn(tables) -> e[a, b]."""

    def __init__(self, c: Circuit, data, tlayers, q, final, e_fn, proxy=False, zgroups=None, gauge_gens=None, xgroups=None, free_tags=()):
        self.c, self.d, self.data = c, c.d, list(data)
        self.n = len(self.data)
        self.pos = {qq: i for i, qq in enumerate(self.data)}
        self.tl = {k: {self.pos[qq]: t for qq, t in v.items()} for k, v in tlayers.items()}
        self.tab = Tables(self.d, q, proxy=proxy)
        ds = set(self.data)
        ffs = [k for k, op in enumerate(c.ops) if op[0] == "FF"]
        rst = [k for k, op in enumerate(c.ops) if op[0] == "R" and k > 0 and set(op[1]) & ds]
        for k in rst:
            assert set(op for op in c.ops[k][1]) <= ds
        self.ff, self.rst = set(ffs), set(rst)
        self.flush = sorted(set(self.tl) | self.ff | self.rst)
        for k, op in enumerate(c.ops):
            if op[0] == "M":
                assert not ds & set(op[1]), "data measured"
            if op[0] == "FF":
                assert all(qq in ds for _, qq, _, _ in op[1]), "FF off data"
        self.bounds = [0] + self.flush + [len(c.ops)]
        self.dgates = []
        for s in range(len(self.bounds) - 1):
            g = []
            for k in range(self.bounds[s], self.bounds[s + 1]):
                op = c.ops[k]
                if k in self.tl or k == 0:
                    continue
                if op[0] in ("F", "FD", "S"):
                    if set(op[1]) <= ds:
                        g.append(op)
                    else:
                        assert not set(op[1]) & ds, "mixed single-qudit op"
                elif op[0] == "CSUM":
                    if all(a in ds and b in ds for a, b in op[1]):
                        g.append(op)
                    else:
                        assert all(not (a in ds and b in ds) for a, b in op[1]), "mixed CSUM layer"
            self.dgates.append(g)
        self.gauge_gens = gauge_gens
        self.free_tags = set(free_tags)
        self._observables()
        self._random()
        self._zred(zgroups)
        self.xgens = [(None if xgroups is None else xgroups(s, self)) for s in range(len(self.bounds) - 1)]
        self.xgens = [None if g is None else [np.asarray(r, np.int64) % self.d for r in g] for g in self.xgens]
        self.final = final
        self.e = e_fn(self)
        self._final_prep()

    def _observables(self):
        n, d, c = self.n, self.d, self.c
        starts = [0] + self.flush
        B = 2 * n * len(starts)
        forced = {}
        for s, k in enumerate(starts):
            for j in range(n):
                for typ in range(2):
                    b = (s * n + j) * 2 + typ
                    mask = np.zeros(B, bool); mask[b] = True
                    P = np.zeros((B, 1, 2), np.int64); P[b, 0, typ] = 1
                    forced.setdefault(k, []).append((mask, [self.data[j]], P))
        fl = frame_pass(c, self.data, set(self.flush), B, np.random.default_rng(0), scale=0.0, forced=forced, skip_ff=True)["flips"]
        fl = fl.reshape(len(starts), n, 2, c.n_meas)
        self.seg_of, self.U, self.V = {}, {}, {}
        for m in range(c.n_meas):
            segs = [s for s in range(len(starts)) if fl[s, :, :, m].any()]
            assert len(segs) <= 1, f"measurement {m} couples in segments {segs}"
            if segs:
                s = segs[0]
                self.seg_of[m] = s
                self.U[m] = (-fl[s, :, 1, m]) % d; self.V[m] = fl[s, :, 0, m] % d
        self.meas_in = [[m for m in self.seg_of if self.seg_of[m] == s] for s in range(len(starts))]

    def _random(self):
        c, d = self.c, self.d
        s = FrameSim(c, 512, np.random.default_rng(7), gauge=True).run(scale=0.0)
        self.rnd = [m for m in range(c.n_meas) if s.meas[:, m].any()]
        det = {}
        self.free = set()                                    # recorded, not post-selected (incoherent split by outcome)
        for j, t in enumerate(c.detectors):
            assert len(t) == 1 and t[0][1] == 1, "multi-term detector"
            if c.det_tags[j] in self.free_tags:
                self.free.add(t[0][0])
            else:
                det[t[0][0]] = j
        for m in range(c.n_meas):
            assert (m in det) + (m in self.rnd) + (m in self.free) == 1, f"measurement {m}: detector / random mismatch"
        for m in self.rnd:
            assert np.gcd.reduce(np.append(s.meas[:, m] % d, d)) == 1, "random outcome not uniform on Z_d"
        self.trivial = [m for m in range(c.n_meas) if m not in self.seg_of and m not in self.rnd and m not in self.free]
        self.E, self.G = {}, {}
        if not self.rnd:
            return
        # FF Paulis per random measurement (data positions)
        assert len(self.ff) == 1
        kff = next(iter(self.ff))
        for m in self.rnd:
            ex = np.zeros(self.n, np.int64); ez = np.zeros(self.n, np.int64)
            self.E[m] = (ex, ez)
        for i, qq, ax, cf in c.ops[kff][1]:
            (self.E[i][0] if ax == "x" else self.E[i][1])[self.pos[qq]] += cf
        for m in self.rnd:
            self.E[m] = (self.E[m][0] % d, self.E[m][1] % d)
        # G_m = Z^{g} on the fresh qudits at the reset segment's start: flips m by 1, other random ones by 0
        segs = {self.seg_of.get(m) for m in self.rnd}
        assert len(segs) == 1 and None not in segs, "random measurements must couple in one segment"
        s = segs.pop()
        kr = self.bounds[s]
        assert kr in self.rst, "random measurements are not in the growth segment"
        self.rseg, self.kff = s, kff
        new = [self.pos[qq] for qq in self.c.ops[kr][1]]
        gens = []                                            # stabilisers of psi at the segment start: Z on fresh qudits,
        for j in new:                                        # plus the old code's stabilisers (gauge_gens)
            gz = np.zeros(self.n, np.int64); gz[j] = 1; gens.append((np.zeros(self.n, np.int64), gz))
        gens += list(self.gauge_gens or [])
        A = np.array([[(self.V[mp] @ gx - self.U[mp] @ gz) % d for gx, gz in gens] for mp in self.rnd])
        for m in self.rnd:
            b = np.array([1 if mp == m else 0 for mp in self.rnd])
            w = _solve_affine(A, b, d)
            assert w is not None, "no gauge stabiliser for a random measurement"
            gx = sum(int(wi) * g_[0] for wi, g_ in zip(w, gens)) % d
            gz = sum(int(wi) * g_[1] for wi, g_ in zip(w, gens)) % d
            for mp in self.meas_in[s]:                      # deterministic ones are not flipped
                if mp not in self.rnd:
                    assert (self.V[mp] @ gx - self.U[mp] @ gz) % d == 0
            self.G[m] = (gx, gz)
        ms = list(self.rnd)
        for i in range(len(ms)):                             # the G_m commute (one abelian stabiliser group)
            for j in range(i):
                (ax, az), (bx, bz) = self.G[ms[i]], self.G[ms[j]]
                assert (az @ bx - bz @ ax) % d == 0
        self.newpos = new

    def _zred(self, zgroups):
        """zgroups(s, plan) -> Z-stabiliser rows (data positions) valid at segment s's start, or None."""
        d = self.d
        self.zred = []
        for s in range(len(self.bounds) - 1):
            M = zgroups(s, self) if zgroups else None
            if M is None or not len(M):
                self.zred.append(None); continue
            W = howell_form(np.asarray(M) % d, d)
            red = []
            for row in W:
                c = int(np.nonzero(row)[0][0]); red.append((c, int(row[c]), np.asarray(row, np.int64)))
            self.zred.append(red)

    def gates(self, rows, ops):
        tab = self.tab
        for op in ops:
            name = op[0]
            if name in ("F", "FD"):
                T = tab.conj(name)
                for qq in op[1]:
                    rows.conj1(self.pos[qq], T)
            elif name == "S":
                for qq, a in zip(op[1], op[2]):
                    rows.conj1(self.pos[qq], tab.conj(("S", int(a) % self.d)))
            elif name == "CSUM":
                for (a, b), w in zip(op[1], op[2]):
                    rows.conj2(self.pos[a], self.pos[b], tab.conj(("CSUM", int(w) % self.d)))

    # --------------------------------------------------------------------------------------- final evaluation
    def _final_prep(self):
        d = self.d
        f = self.final
        xl, zl = np.asarray(f["xl"]) % d, np.asarray(f["zl"]) % d
        H, MZ = f["H"], f["MZ"]
        assert not np.any((MZ @ xl) % d) and not np.any((H @ zl) % d)
        cc = int(zl @ xl) % d
        assert np.gcd(cc, d) == 1, "logical pair not conjugate"
        self.cinv = pow(cc, -1, d); self.cxz = cc
        # every element of ker MZ is rowspan(H) + a xl (and dually): check on kernel generators
        Hr = [list(r) for r in (H % d)]; Mr = [list(r) for r in (MZ % d)]
        for kv in right_kernel_howell(MZ % d, d):
            kv = np.array(kv) % d; a = (int(zl @ kv) * self.cinv) % d
            assert in_span(list((kv - a * xl) % d), Hr, d), "X normaliser not stab + logical"
        for kv in right_kernel_howell(H % d, d):
            kv = np.array(kv) % d; b = (int(xl @ kv) * self.cinv) % d
            assert in_span(list((kv - b * zl) % d), Mr, d), "Z normaliser not stab + logical"
        tx, tz = f["tables"]
        self.tx, self.tz = tx, tz
        e = self.e
        w = np.exp(2j * np.pi / d)
        Gm = np.zeros((d * d, d * d), complex)
        for a in range(d):
            for b in range(d):
                for a2 in range(d):
                    for b2 in range(d):
                        Gm[a * d + b, a2 * d + b2] = w ** ((-b * (a2 - a) * cc) % d) * e[(a2 - a) % d, (b2 - b) % d]
        assert np.allclose(Gm, Gm.conj().T, atol=1e-9) and abs(Gm[0, 0] - 1) < 1e-9
        self.Gram = Gm
        self.xl, self.zl = xl, zl

    def evaluate(self, rows, nshots):
        """-> per shot (w, F_raw, w_triv, F_triv)."""
        d, n = self.d, self.n
        H, MZ = self.final["H"], self.final["MZ"]
        out = np.zeros((nshots, 4))
        if not len(rows):
            return out
        X = rows.X.astype(np.int64); Z = rows.Z.astype(np.int64)
        sz = (X @ MZ.T) % d; sx = (Z @ H.T) % d
        pw_z = d ** np.arange(MZ.shape[0]); pw_x = d ** np.arange(H.shape[0])
        kz = sz @ pw_z; kx = sx @ pw_x
        # sector representative R = X^e Z^e': decoder corrections, or (undecodable) the first row of the sector
        grp_key = np.stack([rows.sid, rows.cls, kz, kx], 1)
        uk, first, inv = np.unique(grp_key, axis=0, return_index=True, return_inverse=True)
        inv = inv.ravel()
        E = np.zeros((len(uk), n), np.int64); Ep = np.zeros((len(uk), n), np.int64); dec = np.zeros(len(uk), bool)
        for g in range(len(uk)):
            i0 = first[g]
            ez = self.tx.get(tuple(sz[i0])) if sz[i0].any() else np.zeros(n, np.int64)
            ex = self.tz.get(tuple(sx[i0])) if sx[i0].any() else np.zeros(n, np.int64)
            if ez is not None and ex is not None:
                E[g], Ep[g], dec[g] = np.asarray(ez) % d, np.asarray(ex) % d, True
            else:
                E[g], Ep[g] = X[i0], Z[i0]
        e_r, ep_r = E[inv], Ep[inv]
        dx = (X - e_r) % d; dz = (Z - ep_r) % d
        ph = np.exp(2j * np.pi * ((-(ep_r * dx).sum(1)) % d) / d)
        a_ = ((dx @ self.zl) * self.cinv) % d; b_ = ((dz @ self.xl) * self.cinv) % d
        beta = np.zeros((len(uk), d * d), complex)
        np.add.at(beta, (inv, a_ * d + b_), rows.a * ph)
        wg = np.einsum("gi,ij,gj->g", beta.conj(), self.Gram, beta).real
        Fg = np.abs(beta @ self.e.reshape(-1)) ** 2
        sid = uk[:, 0]
        triv = (uk[:, 2] == 0) & (uk[:, 3] == 0)
        np.add.at(out[:, 0], sid, wg)
        np.add.at(out[:, 1], sid, Fg * dec)
        np.add.at(out[:, 2], sid, wg * triv)
        np.add.at(out[:, 3], sid, Fg * triv)
        return out

    # --------------------------------------------------------------------------------------- the branch pass
    def run_rows(self, fp, sel, cap=None):
        """fp: frame_pass output; sel: shot indices to simulate.  Returns (eval (len(sel), 4), max rows)."""
        d, n = self.d, self.n
        flips = fp["flips"][sel]
        ok = np.ones(len(sel), bool)
        for m in self.trivial:
            ok &= flips[:, m] == 0
        rows = Rows(n, np.nonzero(ok)[0], d)
        maxr = len(rows)
        for s in range(len(self.bounds) - 1):
            k = self.bounds[s]
            if s > 0:
                if k in self.ff and self.rnd:
                    self._feedforward(rows)
                xr, zr = fp["rec"][k]
                rows.left(xr[sel][rows.sid], zr[sel][rows.sid])
                if k in self.rst:
                    for j in [self.pos[qq] for qq in self.c.ops[k][1]]:
                        assert not rows.X[:, j].any() and not rows.Z[:, j].any(), "branch on a reset qudit"
                if k in self.tl:
                    if self.xgens[s] is not None:
                        rows.reduce_x(self.xgens[s]); rows.merge(self.zred[s])
                    for j, t in self.tl[k].items():
                        rows.texpand(j, self.tab.tcoef(t))
                        rows.merge(self.zred[s])
                        if cap and len(rows) > cap:
                            raise OverflowError(len(rows))
            for m in self.meas_in[s]:
                sh = (rows.X.astype(np.int64) @ self.V[m] - rows.Z.astype(np.int64) @ self.U[m]) % d
                if m in self.rnd:
                    rows.shift[m] = sh.astype(np.int16)
                elif m in self.free:
                    rows.cls = rows.cls * 1000003 + (sh + flips[rows.sid, m]) % d
                else:
                    keep = (sh + flips[rows.sid, m]) % d == 0
                    if not keep.all():
                        rows.take(np.nonzero(keep)[0])
            if s == getattr(self, "rseg", -1):
                for m in self.rnd:                          # P -> P G^{-shift} at the growth segment start
                    sh = rows.shift.get(m)
                    if sh is None:
                        continue
                    gx, gz = self.G[m]
                    k = ((-(sh.astype(np.int64))) % d)
                    # (X^gx Z^gz)^k = w^{(gz.gx) k(k-1)/2} X^{k gx} Z^{k gz}
                    rows.a = rows.a * np.exp(2j * np.pi * (((int(gz @ gx) % d) * (k * (k - 1) // 2)) % d) / d)
                    rows.right((k[:, None] * gx[None, :]) % d, (k[:, None] * gz[None, :]) % d)
            self.gates(rows, self.dgates[s])
            maxr = max(maxr, len(rows))
            if cap and len(rows) > cap:
                raise OverflowError(len(rows))
        xr, zr = fp["rec"]["end"]
        rows.left(xr[sel][rows.sid], zr[sel][rows.sid])
        rows.merge(self.zred[-1])
        return self.evaluate(rows, len(sel)), maxr

    def _feedforward(self, rows):
        """dephasing of the outcome-record sum: class label by s_m = <E_m, P> (P after the growth segment's gates)."""
        d = self.d
        if not len(rows):
            return
        lab = np.zeros(len(rows), np.int64)
        X = rows.X.astype(np.int64); Z = rows.Z.astype(np.int64)
        for m in self.rnd:
            ex, ez = self.E[m]
            s = (ez @ X.T - Z @ ex) % d                     # <E, P>: E P = w^{<E,P>} P E, <A,B> = z_A.x_B - z_B.x_A
            lab = lab * d + s
        rows.cls = rows.cls * (d ** len(self.rnd)) + lab
        rows.shift = {}
        rows.merge()


# ----------------------------------------------------------------------------------------------- e(a, b): <psi|X^{a xl} Z^{b zl}|psi>
def encoder_expectation(plan: Plan, xl, zl, stage_pos=None):
    """backward conjugation of N = X^{a xl} Z^{b zl} through segment-1 gates, the injection T, segment-0 gates, onto |0>.
    stage_pos: data positions of the injected stage (default all)."""
    d, n = plan.d, plan.n
    kinj = plan.flush[0]
    assert kinj in plan.tl and len(plan.tl[kinj]) == 1
    (jinj, tinj), = plan.tl[kinj].items()
    inv = {"F": "FD", "FD": "F"}

    def back(rows, ops):
        for op in reversed(ops):
            name = op[0]
            if name in ("F", "FD"):
                T = plan.tab.conj(inv[name])
                for qq in op[1]:
                    rows.conj1(plan.pos[qq], T)
            elif name == "S":
                for qq, a in zip(op[1], op[2]):
                    rows.conj1(plan.pos[qq], plan.tab.conj(("S", (-int(a)) % d)))
            elif name == "CSUM":
                for (a, b), w in reversed(list(zip(op[1], op[2]))):
                    rows.conj2(plan.pos[a], plan.pos[b], plan.tab.conj(("CSUM", (-int(w)) % d)))
    e = np.zeros((d, d), complex)
    later = plan.dgates[1]
    for a in range(d):
        for b in range(d):
            rows = Rows(n, [0], d)
            rows.X[0] = (a * np.asarray(xl)) % d; rows.Z[0] = (b * np.asarray(zl)) % d
            back(rows, later)
            rows.texpand(jinj, plan.tab.tcoef(-tinj))
            back(rows, plan.dgates[0])
            sel = ~rows.X.any(1)                             # <0| X^x Z^z |0> = delta_{x, 0}
            e[a, b] = rows.a[sel].sum()
    return e


# ----------------------------------------------------------------------------------------------- circuits
CONFIGS = {"par-hooks-f1": (dict(hooks=True, parallel=True, final_block=True), 1.0),
           "par-hooks": (dict(hooks=True, parallel=True, final_block=True), 0.0)}


def _tlayers(c, inj_q, taumap):
    """injection T on q*; each check's S pair -> T^{-tau}, T^{+tau} (signed tau per data qudit, as exact_t_cult)."""
    S_ops = [k for k, op in enumerate(c.ops) if op[0] == "S"]
    tl = {S_ops[0]: {inj_q: 1}}
    assert c.ops[S_ops[0]][1] == [inj_q] and (len(S_ops) - 1) % 2 == 0
    for i in range(1, len(S_ops), 2):
        k1, k2 = S_ops[i], S_ops[i + 1]
        assert c.ops[k1][1] == c.ops[k2][1]
        tl[k1] = {qq: -int(taumap[qq]) for qq in c.ops[k1][1]}
        tl[k2] = {qq: int(taumap[qq]) for qq in c.ops[k2][1]}
    return tl


def _a1_points(tl):
    """first T layer of each check (the ideal state is the code state just before it)."""
    ks = sorted(tl)
    return set(ks[1::2])


def colorcult_plan(d, p, distance=3, cfg="par-hooks-f1", proxy=False, rounds=1, checks=2, free_tags=()):
    kw, f = CONFIGS[cfg]
    cc = ColorCult(distance, d, p, rounds=rounds, checks_per_round=checks, **kw)
    if f:
        cc.c = add_tick_idle(cc.c, f * p)
    tl = _tlayers(cc.c, cc.data[cc.qstar], {cc.data[j]: int(cc.tau[j]) for j in range(cc.n)})
    n = cc.n
    xl = np.ones(n, np.int64); zl = (cc.sigma * cc.tau) % d
    zl = (np.asarray(cc.tau) % d)
    final = dict(H=cc.H % d, MZ=cc.MZ % d, xl=xl, zl=zl, tables=cc._tables())
    MZ = cc.MZ % d
    plan = Plan(cc.c, cc.data, tl, escaping_carrier(d), final, lambda P: encoder_expectation(P, xl, zl), proxy=proxy,
                zgroups=lambda s, P: MZ if s >= 2 else None,
                xgroups=lambda s, P: (cc.H % d) if P.bounds[s] in _a1_points(tl) else None, free_tags=free_tags)
    return cc, plan


def growcult_plan(d, p, cfg="par-hooks-f1", proxy=False, free_tags=()):
    kw, f = CONFIGS[cfg]
    g = GrowCult(d, p, 1, 1, rounds3=1, rounds5=1, checks3=2, checks5=2, **kw)
    if f:
        g.c = add_tick_idle(g.c, f * p)
    emb = list(g.emb)
    d3 = [g.data5[i] for i in emb]
    tl = _tlayers(g.c, d3[g.qstar], {g.data5[j]: int(g.tau5[j]) for j in range(len(g.data5))})
    xbar, zbar = np.asarray(g.xbar) % d, np.asarray(g.zbar) % d
    new = [i for i in range(len(g.data5)) if i not in emb]
    assert all(xbar[i] == 0 for i in g.zset) and all(zbar[i] == 0 for i in new if i not in g.zset)
    final = dict(H=g.H % d, MZ=g.MZ % d, xl=xbar, zl=zbar, tables=g._tables())

    def efn(P):
        xs = np.zeros(P.n, np.int64); zs = np.zeros(P.n, np.int64)
        xs[emb] = xbar[emb]; zs[emb] = zbar[emb]                 # new qudits: X on |+>, Z on |0> -> 1
        return encoder_expectation(P, xs, zs)
    nb = len(g.data5)
    MZ3 = np.zeros((g.MZs.shape[0], nb), np.int64); MZ3[:, emb] = g.MZs % d
    MZ5 = g.MZb % d

    def zg(s, P):
        sff = P.bounds.index(P.kff)
        if 2 <= s < P.rseg:
            return MZ3
        if s >= sff:
            return MZ5
        return None
    gg = []
    for r in g.Hs % d:
        gx = np.zeros(nb, np.int64); gx[emb] = r; gg.append((gx, np.zeros(nb, np.int64)))
    for r in g.MZs % d:
        gz = np.zeros(nb, np.int64); gz[emb] = r; gg.append((np.zeros(nb, np.int64), gz))
    H3 = np.zeros((g.Hs.shape[0], nb), np.int64); H3[:, emb] = g.Hs % d
    a1 = _a1_points(tl)

    def xg(s, P):
        k = P.bounds[s]
        if k not in a1:
            return None
        return H3 if set(P.tl[k]) <= set(emb) else (g.Hb % d)
    plan = Plan(g.c, g.data5, tl, escaping_carrier(d), final, efn, proxy=proxy, zgroups=zg, gauge_gens=gg, xgroups=xg, free_tags=free_tags)
    return g, plan


def g_qstar(g):
    return g.qstar


# ----------------------------------------------------------------------------------------------- Monte Carlo
def proxy_outcome(cc, plan, fp):
    det = np.zeros(fp["measp"].shape[0], bool)
    for j, t in enumerate(cc.c.detectors):
        (m, cf), = t
        if cc.c.det_tags[j] in plan.free_tags:
            continue
        det |= (cf * fp["measp"][:, m]) % plan.d != 0
    okp = ~det
    fail, nontriv = cc.logical_failures(fp["xp"][okp], fp["zp"][okp])
    return okp, fail, nontriv


def run(d, p, shots, out, seed=0, batch=20000, cfg="par-hooks-f1", circuit="grow", sub=200, cap=None, free_tags=()):
    cap = cap or int(os.environ.get("BT_CAP", 400_000))
    if circuit == "grow":
        cc, plan = growcult_plan(d, p, cfg, free_tags=free_tags)
    else:                                                    # "color3" / "color5": ColorCult at that distance
        cc, plan = colorcult_plan(d, p, int(circuit[-1]), cfg, free_tags=free_tags)
    rng = np.random.default_rng(seed)
    tot = dict(d=d, p=p, cfg=cfg, circuit=circuit, seed=seed, shots=0, branch_shots=0, w=0.0, Fraw=0.0, w_triv=0.0, Ftriv=0.0,
               proxy_acc=0, proxy_fail=0, proxy_acc_triv=0, proxy_fail_triv=0, max_rows=0, overflow=0, seconds=0.0,
               inf_on_proxy_ok=0.0, inf_triv_on_proxy_ok=0.0, w_on_proxy_rej=0.0, w_on_proxy_fail=0.0)
    if os.path.exists(out):
        tot = json.load(open(out)); rng = np.random.default_rng([seed, tot["shots"]])
        print("resuming at", tot["shots"], flush=True)
    t0 = time.time() - tot.get("seconds", 0.0)
    while tot["shots"] < shots:
        B = min(batch, shots - tot["shots"])
        fp = frame_pass(cc.c, plan.data, set(plan.flush), B, rng)
        okp, fail, nontriv = proxy_outcome(cc, plan, fp)
        pgood = np.zeros(B, bool); pgood[np.nonzero(okp)[0][~fail]] = True     # proxy: accepted and correct
        pfail = np.zeros(B, bool); pfail[np.nonzero(okp)[0][fail]] = True
        tot["proxy_acc"] += int(okp.sum()); tot["proxy_fail"] += int(fail.sum())
        tot["proxy_acc_triv"] += int((~nontriv).sum()); tot["proxy_fail_triv"] += int((fail & ~nontriv).sum())
        act = np.any(fp["flips"], axis=1)
        for k, (xr, zr) in fp["rec"].items():
            act |= np.any(xr, axis=1) | np.any(zr, axis=1)
        idle = ~act
        tot["w"] += float(idle.sum()); tot["Fraw"] += float(idle.sum()); tot["w_triv"] += float(idle.sum()); tot["Ftriv"] += float(idle.sum())
        ia = np.nonzero(act)[0]
        for i in range(0, len(ia), sub):
            sel = ia[i:i + sub]
            try:
                ev, mr = plan.run_rows(fp, sel, cap=cap)
            except (OverflowError, MemoryError):
                ev = np.zeros((0, 4)); mr = 0
                for s1 in sel:                               # one by one; a single shot above the cap is counted, not run
                    try:
                        e1, m1 = plan.run_rows(fp, np.array([s1]), cap=cap)
                        ev = np.concatenate([ev, e1]); mr = max(mr, m1)
                    except (OverflowError, MemoryError):
                        tot["overflow"] += 1
                        ev = np.concatenate([ev, np.full((1, 4), np.nan)])
            tot["max_rows"] = max(tot["max_rows"], int(mr))
            ovf = np.isnan(ev[:, 0])
            if ovf.any():                                    # overflow shots: counted, excluded (bounded by their number)
                tot["overflow_proxy_acc"] = tot.get("overflow_proxy_acc", 0) + int(okp[sel][ovf].sum())
                ev = np.nan_to_num(ev, nan=0.0)
            if len(ev) == len(sel):
                g_, f_, r_ = pgood[sel], pfail[sel], ~okp[sel]
                tot["inf_on_proxy_ok"] += float((ev[g_, 0] - ev[g_, 1]).sum())
                tot["inf_triv_on_proxy_ok"] += float((ev[g_, 2] - ev[g_, 3]).sum())
                tot["w_on_proxy_rej"] += float(ev[r_, 0].sum()); tot["w_on_proxy_fail"] += float(ev[f_, 0].sum())
            tot["w"] += float(ev[:, 0].sum()); tot["Fraw"] += float(ev[:, 1].sum())
            tot["w_triv"] += float(ev[:, 2].sum()); tot["Ftriv"] += float(ev[:, 3].sum())
        tot["shots"] += B; tot["branch_shots"] += int(len(ia))
        tot["seconds"] = time.time() - t0
        tot["accept_exact"] = tot["w"] / tot["shots"]; tot["eps_exact"] = 1 - tot["Fraw"] / tot["w"]
        tot["accept_proxy"] = tot["proxy_acc"] / tot["shots"]; tot["eps_proxy"] = tot["proxy_fail"] / max(tot["proxy_acc"], 1)
        tot["eps_exact_triv"] = 1 - tot["Ftriv"] / max(tot["w_triv"], 1e-300)
        tot["eps_proxy_triv"] = tot["proxy_fail_triv"] / max(tot["proxy_acc_triv"], 1)
        json.dump(tot, open(out + ".tmp", "w")); os.replace(out + ".tmp", out)
        print(f"progress {tot['shots']}/{shots} ({tot['seconds'] / 60:.1f} min): accept exact {tot['accept_exact']:.4f} proxy "
              f"{tot['accept_proxy']:.4f}; eps exact {tot['eps_exact']:.3e} proxy {tot['eps_proxy']:.3e} ({tot['proxy_fail']} proxy fails); "
              f"triv exact {tot['eps_exact_triv']:.3e} proxy {tot['eps_proxy_triv']:.3e} ({tot['proxy_fail_triv']}); "
              f"max rows {tot['max_rows']}, overflow {tot['overflow']}", flush=True)
    return tot


if __name__ == "__main__":
    if sys.argv[1] == "run":
        a = sys.argv[2:]
        run(int(a[0]), float(a[1]), int(a[2]), a[3], int(a[4]) if len(a) > 4 else 0, int(a[5]) if len(a) > 5 else 20000,
            a[6] if len(a) > 6 else "par-hooks-f1", a[7] if len(a) > 7 else "grow",
            free_tags=tuple(a[8].split(",")) if len(a) > 8 and a[8] != "-" else ())
