"""Exact-T validation of the Clifford proxy for stage-1 cultivation (ColorCult, D = 3).

Hybrid simulator, exact for Pauli noise.  Every gate of the circuit is Clifford except the T layers on the data qudits.
ColorCult's proxy S^sigma ops (injection; A1 = T^{-tau}, A2 = T^{tau} around each check's coupling) are replaced by the
exact escaping T (cultivation.escaping_carrier).  Data qudits are never measured, and every measurement is a
post-selected detector.  Per sampled Pauli fault set:
  frame pass  : the Z_d Pauli frame of every fault, propagated through the Clifford ops; the data frame is flushed (recorded
                and zeroed) just before every T layer and at the end; measurement flips are recorded.
  dense pass  : the n data qudits as a state vector, segment by segment.  At each T layer: apply the flushed frame, the
                exact T layer, the ideal projective measurement of every detector coupled in that segment, and the
                segment's data-only Clifford gates (encoder).  Each detector is post-selected to its reference: index
                k0 - flip.  Its data observable O (at the segment start) comes from unit X_j / Z_j injections into the
                noiseless frame pass (u_j = -flip(Z_j), v_j = flip(X_j)).
  output      : acceptance weight w = |psi|^2; infidelity after ideal min-weight decoding (ColorCult's lookup tables, as
                the proxy's logical_failures); and the same with final-syndrome post-selection.
The proxy outcome of every shot comes from the same noise draws (paired comparison).

Validation: `selftest` (dense kernels vs frame conjugation rules, noiseless determinism at d = 2, 4, 8, brute-force
full-register dense on toy circuits at d = 4, 8 and on ColorCult D = 3 at d = 2).

    PYTHONPATH=. python scratchpad/exact_t_cult.py selftest
    PYTHONPATH=. python scratchpad/exact_t_cult.py run d p shots out.json [seed] [batch]
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

from primepower_msd.cult_circuit import ColorCult
from primepower_msd.cultivation import escaping_carrier
from primepower_msd.zdsim import Circuit


# ----------------------------------------------------------------------------------------------- dense kernels
class Dense:
    """Batch of states on n qudits of dimension d, array (B,) + (d,)*n, axis 1+j = qudit j."""

    def __init__(self, n, d):
        self.n, self.d = n, d
        self.w = np.exp(2j * np.pi / d)
        self.Fm = np.array([[self.w ** (j * k) for j in range(d)] for k in range(d)]) / np.sqrt(d)   # F[k, j]

    def shape(self, B):
        return (B,) + (self.d,) * self.n

    def X(self, psi, j, a):
        """X^a on qudit j; a scalar or per-shot array."""
        if np.isscalar(a):
            return np.roll(psi, int(a) % self.d, axis=1 + j)
        out = psi.copy()
        a = np.asarray(a) % self.d
        for v in range(1, self.d):
            m = a == v
            if m.any():
                out[m] = np.roll(psi[m], v, axis=1 + j)
        return out

    def Z(self, psi, j, z):
        d = self.d
        z = np.broadcast_to(np.asarray(z) % d, (psi.shape[0],))
        ph = self.w ** (np.outer(z, np.arange(d)) % d)
        sh = [psi.shape[0]] + [1] * self.n; sh[1 + j] = d
        return psi * ph.reshape(sh)

    def diag(self, psi, j, ph):
        sh = [1] * (self.n + 1); sh[1 + j] = self.d
        return psi * np.asarray(ph).reshape(sh)

    def F(self, psi, j, inv=False):
        M = self.Fm.conj().T if inv else self.Fm
        return np.moveaxis(np.tensordot(M, psi, axes=([1], [1 + j])), 0, 1 + j)

    def S(self, psi, j, a):
        d = self.d
        x = np.arange(d)
        return self.diag(psi, j, np.exp(2j * np.pi * (a * x * x % (2 * d)) / (2 * d)))

    def T(self, psi, j, t, q):
        M = 4 * self.d
        return self.diag(psi, j, np.exp(2j * np.pi * ((t * np.asarray(q)) % M) / M))

    def CSUM(self, psi, c, t, w):
        """|a, b> -> |a, b + w a>."""
        out = np.empty_like(psi)
        for a in range(self.d):
            idx = [slice(None)] * psi.ndim; idx[1 + c] = a
            out[tuple(idx)] = np.roll(psi[tuple(idx)], (w * a) % self.d, axis=(1 + t) - (1 if t > c else 0))
        return out

    def pauli(self, psi, x, z):
        """apply X^x Z^z (per-shot (B, n) arrays): Z first, then X."""
        for j in range(self.n):
            if np.any(z[:, j]):
                psi = self.Z(psi, j, z[:, j])
            if np.any(x[:, j]):
                psi = self.X(psi, j, x[:, j])
        return psi


# ----------------------------------------------------------------------------------------------- frame pass
def frame_pass(c: Circuit, data, flush_ops, B, rng, scale=1.0, forced=None):
    """Two frames from the same noise draws: the proxy frame (plain FrameSim semantics, gauge off) and the hybrid frame
    (data part flushed before every op in flush_ops and at the end).  Returns dict with rec[k] = (x, z) data frames,
    flips (B, n_meas), proxy (x, z) final frames and proxy measurements."""
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
            elif name == "FF":
                raise NotImplementedError("FF")
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


# ----------------------------------------------------------------------------------------------- hybrid plan
class Hybrid:
    """c: Circuit; data: data qudit ids; tlayers: {op index: {qudit: t}} (exact T^t replaces that op); q: T table."""

    def __init__(self, c: Circuit, data, tlayers, q):
        self.c, self.d, self.data, self.q = c, c.d, list(data), np.asarray(q)
        self.n = len(self.data)
        self.pos = {qq: i for i, qq in enumerate(self.data)}
        self.tl = {k: {self.pos[qq]: t for qq, t in v.items()} for k, v in tlayers.items()}
        self.flush = sorted(tlayers)
        self.D = Dense(self.n, self.d)
        ds = set(self.data)
        for k, op in enumerate(c.ops):
            if op[0] == "M":
                assert not ds & set(op[1]), "data measured"
            if op[0] == "R" and ds & set(op[1]):
                assert k == 0 and set(op[1]) == ds, "data reset after start"
        # segments: 0 = [0, flush[0]), s = [flush[s-1], flush[s]) starting with its T layer
        self.bounds = [0] + self.flush + [len(c.ops)]
        self.dgates = []          # per segment: data-only Clifford ops (positions)
        for s in range(len(self.bounds) - 1):
            g = []
            for k in range(self.bounds[s], self.bounds[s + 1]):
                op = c.ops[k]
                if k in self.tl:
                    continue
                if op[0] in ("F", "FD", "S"):
                    qs = op[1]
                    if set(qs) <= ds:
                        g.append((k, op))
                    else:
                        assert not set(qs) & ds, "mixed single-qudit op"
                elif op[0] == "CSUM":
                    if all(a in ds and b in ds for a, b in op[1]):
                        g.append((k, op))
                    else:
                        assert all(not (a in ds and b in ds) for a, b in op[1]), "mixed CSUM layer"
            self.dgates.append(g)
        self._observables()
        self._reference()

    # O_m (data Pauli at its segment start) from unit injections
    def _observables(self):
        n, d, c = self.n, self.d, self.c
        starts = [0] + self.flush            # inject after this op (seg 0: after the data reset)
        B = 2 * n * len(starts)
        forced = {}
        for s, k in enumerate(starts):
            for j in range(n):
                for typ in range(2):
                    b = (s * n + j) * 2 + typ
                    mask = np.zeros(B, bool); mask[b] = True
                    P = np.zeros((B, 1, 2), np.int64); P[b, 0, typ] = 1
                    forced.setdefault(k, []).append((mask, [self.data[j]], P))
        fl = frame_pass(c, self.data, set(self.flush), B, np.random.default_rng(0), scale=0.0, forced=forced)["flips"]
        fl = fl.reshape(len(starts), n, 2, c.n_meas)
        self.seg_of, self.O = {}, {}
        for m in range(c.n_meas):
            segs = [s for s in range(len(starts)) if fl[s, :, :, m].any()]
            assert len(segs) <= 1, f"measurement {m} couples in segments {segs}"
            if segs:
                s = segs[0]
                u = (-fl[s, :, 1, m]) % d; v = fl[s, :, 0, m] % d
                self.seg_of[m], self.O[m] = s, (u, v)
        self.meas_in = [[m for m in self.seg_of if self.seg_of[m] == s] for s in range(len(starts))]
        self.trivial = [m for m in range(c.n_meas) if m not in self.seg_of]
        self._ops = {m: self._op_arrays(*self.O[m]) for m in self.O}
        d, n = self.d, self.n
        self.grid = np.indices((d,) * n).reshape(n, -1)
        self.lab = {}                        # m -> ('Z' | 'X', label over comp / Fourier basis); pure observables only
        for m, (u, v) in self.O.items():
            if not u.any():
                self.lab[m] = ("Z", (v @ self.grid) % d)
            elif not v.any():
                self.lab[m] = ("X", (-(u @ self.grid)) % d)
        self.tphase = {k: np.prod([np.exp(2j * np.pi * ((t * self.q[self.grid[j]]) % (4 * d)) / (4 * d)) for j, t in tl.items()], axis=0)
                       for k, tl in self.tl.items()}

    def _op_arrays(self, u, v):
        """O = X^u Z^v as (phase over basis, gather index); normalised so that O^d = 1."""
        d, n = self.d, self.n
        grid = np.indices((d,) * n).reshape(n, -1)
        ph = self.D.w ** ((v @ grid) % d)
        src = np.ravel_multi_index(tuple((grid - u[:, None]) % d), (d,) * n)
        vec = np.random.default_rng(1).normal(size=d ** n) + 0j
        y = vec.copy()
        for _ in range(d):
            y = (ph * y)[src]
        cst = y[0] / vec[0]
        assert np.allclose(y, cst * vec)
        return ph, src, cst ** (-1.0 / d)

    def _apply_O(self, psi_flat, m):
        ph, src, nrm = self._ops[m]
        return nrm * (psi_flat * ph)[:, src]

    def _project(self, psi, m, k):
        """onto eigen-index k (per-shot) of O'_m."""
        B = psi.shape[0]
        y = psi.reshape(B, -1)
        acc = y.copy() / self.d
        cur = y
        for j in range(1, self.d):
            cur = self._apply_O(cur, m)
            acc = acc + (self.D.w ** ((-j * k) % self.d))[:, None] * cur / self.d
        return acc.reshape(psi.shape)

    def _apply_gate(self, psi, op):
        D, name = self.D, op[0]
        if name in ("F", "FD"):
            for qq in op[1]:
                psi = D.F(psi, self.pos[qq], inv=(name == "FD"))
        elif name == "S":
            for qq, a in zip(op[1], op[2]):
                psi = D.S(psi, self.pos[qq], a)
        elif name == "CSUM":
            for (a, b), w in zip(op[1], op[2]):
                psi = D.CSUM(psi, self.pos[a], self.pos[b], w)
        return psi

    def dense(self, rec, flips, k0=None):
        """rec: frame records (from frame_pass, already sub-selected); returns psi (B, d^n) unnormalised."""
        B = flips.shape[0]
        psi = np.zeros(self.D.shape(B), complex)
        psi[(slice(None),) + (0,) * self.n] = 1.0
        ok = np.ones(B, bool)
        for m in self.trivial:
            ok &= flips[:, m] == 0
        for s in range(len(self.bounds) - 1):
            if s > 0:
                k = self.flush[s - 1]
                psi = self.pauli_flat(psi.reshape(B, -1), *rec[k]) * self.tphase[k][None, :]
                psi = psi.reshape(self.D.shape(B))
            if k0 is None:
                for m in self.meas_in[s]:  # reference pass: record the noiseless eigen-index
                    y = psi.reshape(B, -1)
                    ev = np.einsum("bi,bi->b", y.conj(), self._apply_O(y, m))
                    assert np.allclose(np.abs(ev), 1, atol=1e-9), f"meas {m}: noiseless not an eigenstate {ev[:3]}"
                    self.k0[m] = int(np.round(np.angle(ev[0]) / (2 * np.pi / self.d))) % self.d
            else:
                ms = self.meas_in[s]
                y = psi.reshape(B, -1)
                for m in ms:
                    if m in self.lab and self.lab[m][0] == "Z":
                        y = y * (self.lab[m][1][None, :] == ((k0[m] - flips[:, m]) % self.d)[:, None])
                xs = [m for m in ms if m in self.lab and self.lab[m][0] == "X"]
                if xs:
                    sh = self.D.shape(B); ax = tuple(range(1, self.n + 1))
                    yf = np.fft.fftn(y.reshape(sh), axes=ax, norm="ortho").reshape(B, -1)
                    for m in xs:
                        yf = yf * (self.lab[m][1][None, :] == ((k0[m] - flips[:, m]) % self.d)[:, None])
                    y = np.fft.ifftn(yf.reshape(sh), axes=ax, norm="ortho").reshape(B, -1)
                psi = y.reshape(self.D.shape(B))
                for m in ms:
                    if m not in self.lab:
                        psi = self._project(psi, m, (k0[m] - flips[:, m]) % self.d)
            for _, op in self.dgates[s]:
                psi = self._apply_gate(psi, op)
        psi = self.pauli_flat(psi.reshape(B, -1), *rec["end"])
        return psi * ok[:, None]

    def pauli_flat(self, y, x, z):
        """X^x Z^z per shot on flat states (B, d^n)."""
        d = self.d
        if np.any(z):
            y = y * self.D.w ** ((z @ self.grid) % d)
        if np.any(x):
            src = np.ravel_multi_index(tuple((self.grid[:, None, :] - x.T[:, :, None]) % d), (d,) * self.n)   # (B, d^n)
            y = np.take_along_axis(y, src, axis=1)
        return y

    def _reference(self):
        self.k0 = {}
        z = np.zeros((1, self.n), np.int64)
        rec = {k: (z, z) for k in self.flush}; rec["end"] = (z, z)
        self.dense(rec, np.zeros((1, self.c.n_meas), np.int64), k0=None)
        psi = self.dense(rec, np.zeros((1, self.c.n_meas), np.int64), k0=self.k0)
        assert np.isclose(np.vdot(psi[0], psi[0]).real, 1.0), "noiseless run not accepted with certainty"
        self.out0 = psi[0]


# ----------------------------------------------------------------------------------------------- ColorCult wrapper
CONFIGS = {   # the D = 3 stage of cult_escape_e2e's configs (one stabiliser round + 2 checks, as GrowCult's D=3 part), + final round
    "serial": (dict(rounds=2, checks_per_round=2), 0.0),        # CULTIVATION_PLAN section 11 circuit
    "par-hooks": (dict(rounds=1, checks_per_round=2, hooks=True, parallel=True, final_block=True), 0.0),
    "par-hooks-f1": (dict(rounds=1, checks_per_round=2, hooks=True, parallel=True, final_block=True), 1.0),
    "par-hooks-f1-vcat": (dict(rounds=1, checks_per_round=2, hooks=True, parallel=True, final_block=True, verify_cat=True), 1.0),
}


def colorcult_hybrid(d, p, distance=3, cfg="serial", **kw):
    from primepower_msd.zdsim import add_tick_idle
    ckw, f = CONFIGS[cfg]
    cc = ColorCult(distance, d, p, **{**ckw, **kw})
    if f:
        cc.c = add_tick_idle(cc.c, f * p, no_anchor=cc.vanc or ())   # verified cat stays just-in-time
    assert cc.sigma == 1 and cc.inj == 1
    S_ops = [k for k, op in enumerate(cc.c.ops) if op[0] == "S"]
    tl = {S_ops[0]: {cc.data[cc.qstar]: 1}}
    assert cc.c.ops[S_ops[0]][1] == [cc.data[cc.qstar]]
    for i in range(1, len(S_ops), 2):           # check: A1 = T^{-tau}, A2 = T^{tau}
        tl[S_ops[i]] = {cc.data[j]: -int(cc.tau[j]) for j in range(cc.n)}
        tl[S_ops[i + 1]] = {cc.data[j]: int(cc.tau[j]) for j in range(cc.n)}
    return cc, Hybrid(cc.c, cc.data, tl, escaping_carrier(d))


class Decoder:
    """Ideal min-weight decode of ColorCult's code on the dense data register (lookup tables of ColorCult._tables)."""

    def __init__(self, cc: ColorCult, target):
        d, n = cc.d, cc.n
        self.cc, self.d, self.n = cc, d, n
        tx, tz = cc._tables()
        grid = np.indices((d,) * n).reshape(n, -1)
        nz = cc.MZ.shape[0]; nx = cc.H.shape[0]
        self.lab_z = ((cc.MZ @ grid) % d).T @ (d ** np.arange(nz))          # Z-stab syndrome of comp-basis states
        self.lab_x = ((cc.H @ grid) % d).T @ (d ** np.arange(nx))           # X-stab syndrome of Fourier-basis states
        self.D = Dense(n, d)
        T = target.reshape((1,) + (d,) * n)
        self.T = target
        Tf = np.fft.fftn(T, axes=tuple(range(1, n + 1)), norm="ortho")   # coefficients <~k|T>, ~k = F|k>... (see _fourier)
        Tf = self._fourier(T)
        assert np.allclose(np.abs(T.reshape(-1))[self.lab_z != 0], 0) and np.allclose(np.abs(Tf.reshape(-1))[self.lab_x != 0], 0)
        self.sz_list, self.cx = [], []
        for s in range(d ** nz):
            key = tuple((s // d ** r) % d for r in range(nz))
            e = np.zeros(n, np.int64) if s == 0 else tx.get(key)
            if e is not None:
                self.sz_list.append(s); self.cx.append(np.asarray(e))
        V = []
        self.sx_list = []
        for s in range(d ** nx):
            key = tuple((s // d ** r) % d for r in range(nx))
            e = np.zeros(n, np.int64) if s == 0 else tz.get(key)
            if e is not None:
                self.sx_list.append(s)
                y = Tf
                for j in range(n):                 # Z^{e} T in Fourier coordinates = shift k -> k + e
                    y = np.roll(y, int(e[j]) % d, axis=1 + j)
                V.append(y.reshape(-1))
        self.V = np.array(V)

    def _fourier(self, psi):
        """coefficients on |~k> = F|k> = sum_j w^{jk}|j>/sqrt d: c_k = sum_j w^{-jk} psi_j / sqrt d (numpy fft sign)."""
        return np.fft.fftn(psi, axes=tuple(range(1, self.n + 1)), norm="ortho")

    def fidelity(self, psi):
        """psi (B, d^n) unnormalised -> (F_raw, w_triv, F_triv_raw): sum over decodable syndromes of |<T|C_s Pi_s psi>|^2."""
        B = psi.shape[0]
        d, n = self.d, self.n
        Fr = np.zeros(B)
        Vc = self.V.conj().T
        for s, e in zip(self.sz_list, self.cx):
            ph = psi * (self.lab_z == s)[None, :]
            if not np.any(ph):
                continue
            y = ph.reshape((B,) + (d,) * n)
            for j in range(n):
                if e[j]:
                    y = np.roll(y, (-int(e[j])) % d, axis=1 + j)       # X^{-e}
            A = self._fourier(y).reshape(B, -1) @ Vc
            Fr += (np.abs(A) ** 2).sum(1)
        m0 = (self.lab_z == 0)
        y0 = psi * m0[None, :]
        y0f = self._fourier(y0.reshape((B,) + (d,) * n)).reshape(B, -1) * (self.lab_x == 0)[None, :]
        w_triv = (np.abs(y0f) ** 2).sum(1)
        F_triv = np.abs(psi @ self.T.conj()) ** 2
        return Fr, w_triv, F_triv


# ----------------------------------------------------------------------------------------------- brute force (validation)
def brute(c: Circuit, tlayers, q, forced_list, data):
    """full-register dense: every measurement projected on outcome 0; returns the data state (ancillas sliced at 0)."""
    D = Dense(c.nq, c.d)
    psi = np.zeros(D.shape(1), complex); psi[(0,) + (0,) * c.nq] = 1
    for k, op in enumerate(c.ops):
        name = op[0]
        if k in tlayers:
            for qq, t in tlayers[k].items():
                psi = D.T(psi, qq, t, q)
        elif name in ("F", "FD"):
            for qq in op[1]:
                psi = D.F(psi, qq, inv=(name == "FD"))
        elif name == "S":
            for qq, a in zip(op[1], op[2]):
                psi = D.S(psi, qq, a)
        elif name == "CSUM":
            for (a, b), w in zip(op[1], op[2]):
                psi = D.CSUM(psi, a, b, w)
        elif name == "M":
            for qq in op[1]:
                idx = [slice(None)] * psi.ndim; idx[1 + qq] = slice(1, None)
                psi[tuple(idx)] = 0
        elif name == "R":
            if k > 0:
                for qq in op[1]:          # valid for a qudit in a computational basis state, unentangled (post-selected
                    sl = [np.take(psi, v, axis=1 + qq) for v in range(c.d)]           # ancilla, possibly hit by a Pauli)
                    nz = [v for v in range(c.d) if np.any(np.abs(sl[v]) > 1e-12)]
                    assert len(nz) <= 1, "reset of an entangled / superposed qudit"
                    new = np.zeros_like(psi)
                    idx = [slice(None)] * psi.ndim; idx[1 + qq] = 0
                    if nz:
                        new[tuple(idx)] = sl[nz[0]]
                    psi = new
        for (kk, qq, P) in forced_list:
            if kk == k:
                psi = D.Z(psi, qq, P[1]); psi = D.X(psi, qq, P[0])
    for qq in range(c.nq):                # measured ancillas may carry a late Pauli: collapse each onto |0> (basis state)
        if qq in data:
            continue
        sl = [np.take(psi, v, axis=1 + qq) for v in range(c.d)]
        nz = [v for v in range(c.d) if np.any(np.abs(sl[v]) > 1e-12)]
        assert len(nz) <= 1, "final ancilla entangled"
        new = np.zeros_like(psi); idx = [slice(None)] * psi.ndim; idx[1 + qq] = 0
        if nz:
            new[tuple(idx)] = sl[nz[0]]
        psi = new
    idx = [0] + [slice(None) if qq in data else 0 for qq in range(c.nq)]
    sub = psi[tuple(idx)]
    sub = np.transpose(sub, [sorted(data).index(qq) for qq in data])
    return sub.reshape(-1)


def hybrid_forced(H: Hybrid, forced_list):
    B = 1
    forced = {}
    for (k, qq, P) in forced_list:
        forced.setdefault(k, []).append((np.ones(1, bool), [qq], np.array(P).reshape(1, 1, 2)))
    fp = frame_pass(H.c, H.data, set(H.flush), B, np.random.default_rng(0), scale=0.0, forced=forced)
    return H.dense(fp["rec"], fp["flips"], H.k0)[0]


def toy(d):
    """2 data qudits (rep code |T>|.> -> sum_j t_j|jj>), Z-stab, two T-checks with a 2-qudit cat (A1 = T^-1 on q0)."""
    c = Circuit(d)
    q0, q1 = c.q(2); a = c.q(1)[0]; c0, c1 = c.q(2)
    tl = {}
    c.R([q0, q1]); c.F([q0]); c.S([q0], 1); tl[len(c.ops) - 1] = {q0: 1}
    c.CSUM([(q0, q1)], 1)
    for _ in range(2):
        c.R([a]); c.CSUM([(q0, a)], 1); c.CSUM([(q1, a)], -1 % d); m = c.M([a]); c.detector([(m[0], 1)])
        c.R([c0, c1]); c.F([c0]); c.CSUM([(c0, c1)], 1)
        c.S([q0], d - 1); tl[len(c.ops) - 1] = {q0: -1}
        c.CSUM([(c0, q0), (c1, q1)], 1)
        c.S([q0], 1); tl[len(c.ops) - 1] = {q0: 1}
        c.CSUM([(c0, c1)], -1 % d); c.FD([c0]); m = c.M([c0, c1])
        c.F([q1]); c.FD([q1])                       # data-only gates inside a segment
    return c, [q0, q1], tl


def selftest():
    rng = np.random.default_rng(5)
    # 1. dense kernels vs frame conjugation rules: G P G^dag = P' (up to phase)
    for d in (2, 4, 8):
        D = Dense(2, d)
        psi = rng.normal(size=(1, d, d)) + 1j * rng.normal(size=(1, d, d))
        for gate in ("F", "FD", "S", "CSUM"):
            for _ in range(5):
                x, z = rng.integers(0, d, 2), rng.integers(0, d, 2)
                a = int(rng.integers(1, d)); w = int(rng.integers(1, d))
                X1 = np.array([x]); Z1 = np.array([z])
                lhs = D.pauli(psi, X1, Z1)
                if gate == "F": lhs = D.F(lhs, 0); g = lambda s: D.F(s, 0); x2, z2 = x.copy(), z.copy(); x2[0], z2[0] = (-z[0]) % d, x[0]
                elif gate == "FD": lhs = D.F(lhs, 0, inv=True); g = lambda s: D.F(s, 0, inv=True); x2, z2 = x.copy(), z.copy(); x2[0], z2[0] = z[0], (-x[0]) % d
                elif gate == "S": lhs = D.S(lhs, 0, a); g = lambda s: D.S(s, 0, a); x2, z2 = x.copy(), z.copy(); z2[0] = (z[0] + a * x[0]) % d
                else: lhs = D.CSUM(lhs, 0, 1, w); g = lambda s: D.CSUM(s, 0, 1, w); x2, z2 = x.copy(), z.copy(); x2[1] = (x[1] + w * x[0]) % d; z2[0] = (z[0] - w * z[1]) % d
                rhs = D.pauli(g(psi), np.array([x2]), np.array([z2]))
                ov = abs(np.vdot(lhs.ravel(), rhs.ravel())) / (np.linalg.norm(lhs) * np.linalg.norm(rhs))
                assert np.isclose(ov, 1), (d, gate, x, z, ov)
    print("kernels vs frame rules: ok", flush=True)
    # 2. toy circuits, brute force vs hybrid, random forced faults
    for d in (4, 8):
        c, data, tl = toy(d)
        H = Hybrid(c, data, tl, escaping_carrier(d))
        print(f"toy d={d}: segments {len(H.bounds)-1}, O per meas {len(H.O)}, trivial {len(H.trivial)}", flush=True)
        worst = 0
        for trial in range(300):
            nf = int(rng.integers(1, 4))
            fl = [(int(rng.integers(0, len(c.ops))), int(rng.integers(0, c.nq)), (int(rng.integers(0, d)), int(rng.integers(0, d))))
                  for _ in range(nf)]
            b = brute(c, tl, escaping_carrier(d), fl, data)
            h = hybrid_forced(H, fl)
            nb, nh = np.vdot(b, b).real, np.vdot(h, h).real
            ov = abs(np.vdot(b, h)) ** 2
            err = abs(nb - nh) + abs(ov - nb * nh)
            worst = max(worst, err)
            assert err < 1e-9, (d, fl, nb, nh, ov)
        print(f"toy d={d}: brute == hybrid on 300 random fault sets (max dev {worst:.1e})", flush=True)
    # 3. ColorCult noiseless determinism and decode of the noiseless output
    for d in (2, 4, 8):
        if d == 8:
            print("ColorCult d=8 skipped in selftest (8^7 register)", flush=True); continue
        cc, H = colorcult_hybrid(d, 0.0)
        dec = Decoder(cc, H.out0)
        Fr, wt, Ft = dec.fidelity(H.out0[None, :])
        assert np.isclose(Fr[0], 1) and np.isclose(Ft[0], 1)
        print(f"ColorCult D=3 d={d}: noiseless accepted with certainty, decoded fidelity {Fr[0]:.12f}; "
              f"meas with data observable {len(H.O)}, trivial {len(H.trivial)}", flush=True)
        # single data errors after the circuit are corrected
        for j in range(cc.n):
            for P in ((1, 0), (0, 1), (1, 1)):
                rec0 = {k: (np.zeros((1, cc.n), np.int64),) * 2 for k in H.flush}
                xe = np.zeros((1, cc.n), np.int64); ze = np.zeros_like(xe); xe[0, j], ze[0, j] = P
                rec0["end"] = (xe, ze)
                psi = H.dense(rec0, np.zeros((1, cc.c.n_meas), np.int64), H.k0)
                assert np.isclose(dec.fidelity(psi)[0][0], 1), (d, j, P)
        print(f"  single final data errors decoded with fidelity 1", flush=True)
    # 4. ColorCult D=3 d=2: brute force on all 20 qubits vs hybrid
    cc, H = colorcult_hybrid(2, 0.0)
    worst = 0
    for trial in range(40):
        nf = int(rng.integers(1, 4))
        fl = [(int(rng.integers(0, len(cc.c.ops))), int(rng.integers(0, cc.c.nq)), (int(rng.integers(0, 2)), int(rng.integers(0, 2))))
              for _ in range(nf)]
        b = brute(cc.c, H.tl_ids if hasattr(H, "tl_ids") else {k: {H.data[j]: t for j, t in v.items()} for k, v in H.tl.items()},
                  escaping_carrier(2), fl, cc.data)
        h = hybrid_forced(H, fl)
        nb, nh = np.vdot(b, b).real, np.vdot(h, h).real
        ov = abs(np.vdot(b, h)) ** 2
        err = abs(nb - nh) + abs(ov - nb * nh)
        worst = max(worst, err)
        assert err < 1e-9, (fl, nb, nh, ov)
    print(f"ColorCult D=3 d=2: brute (2^20) == hybrid on 40 random fault sets (max dev {worst:.1e})", flush=True)
    print("selftest passed", flush=True)


# ----------------------------------------------------------------------------------------------- Monte Carlo
def run(d, p, shots, out, seed=0, batch=20000, sub=128, cfg="serial"):
    cc, H = colorcult_hybrid(d, p, cfg=cfg)
    dec = Decoder(cc, H.out0)
    rng = np.random.default_rng(seed)
    tot = dict(d=d, p=p, distance=3, cfg=cfg, seed=seed, shots=0, dense_shots=0, w=0.0, Fraw=0.0, w_triv=0.0, Ftriv=0.0,
               proxy_acc=0, proxy_fail=0, proxy_acc_triv=0, proxy_fail_triv=0, w_on_proxy_fail=0.0, inf_on_proxy_ok=0.0)
    if os.path.exists(out):
        tot = json.load(open(out)); rng = np.random.default_rng([seed, tot["shots"]])
        print("resuming at", tot["shots"], flush=True)
    t0 = time.time()
    while tot["shots"] < shots:
        B = min(batch, shots - tot["shots"])
        fp = frame_pass(cc.c, cc.data, set(H.flush), B, rng)
        # proxy outcome (same noise)
        okp = ~np.any(fp["measp"], axis=1)
        fail_p, nontriv = cc.logical_failures(fp["xp"][okp], fp["zp"][okp])
        tot["proxy_acc"] += int(okp.sum()); tot["proxy_fail"] += int(fail_p.sum())
        tot["proxy_acc_triv"] += int((~nontriv).sum()); tot["proxy_fail_triv"] += int((fail_p & ~nontriv).sum())
        pf = np.zeros(B, bool); pf[np.nonzero(okp)[0][fail_p]] = True
        # exact: shots with no recorded frame / flip are noiseless
        act = np.any(fp["flips"], axis=1)
        for k, (xr, zr) in fp["rec"].items():
            act |= np.any(xr, axis=1) | np.any(zr, axis=1)
        idle = ~act
        for m in H.trivial:                   # a flipped flag-type detector rejects the shot with certainty (weight 0)
            act &= fp["flips"][:, m] == 0
        tot["w"] += float(idle.sum()); tot["Fraw"] += float(idle.sum())
        tot["w_triv"] += float(idle.sum()); tot["Ftriv"] += float(idle.sum())
        ia = np.nonzero(act)[0]
        for i in range(0, len(ia), sub):
            sel = ia[i:i + sub]
            rec = {k: (xr[sel], zr[sel]) for k, (xr, zr) in fp["rec"].items()}
            psi = H.dense(rec, fp["flips"][sel], H.k0)
            w = (np.abs(psi) ** 2).sum(1)
            live = w > 1e-14
            if not live.any():
                continue
            Fr, wt, Ft = dec.fidelity(psi[live])
            tot["w"] += float(w.sum()); tot["Fraw"] += float(Fr.sum())
            tot["w_triv"] += float(wt.sum()); tot["Ftriv"] += float(Ft.sum())
            wl = w[live]; s_l = sel[live]
            tot["w_on_proxy_fail"] += float(wl[pf[s_l]].sum())
            tot["inf_on_proxy_ok"] += float((wl - Fr)[~pf[s_l]].sum())
        tot["shots"] += B; tot["dense_shots"] += int(len(ia))
        tot["accept_exact"] = tot["w"] / tot["shots"]; tot["eps_exact"] = 1 - tot["Fraw"] / tot["w"]
        tot["accept_proxy"] = tot["proxy_acc"] / tot["shots"]; tot["eps_proxy"] = tot["proxy_fail"] / max(tot["proxy_acc"], 1)
        tot["eps_exact_triv"] = 1 - tot["Ftriv"] / max(tot["w_triv"], 1e-300)
        tot["eps_proxy_triv"] = tot["proxy_fail_triv"] / max(tot["proxy_acc_triv"], 1)
        tot["seconds"] = tot.get("seconds", 0) + 0
        json.dump(tot, open(out + ".tmp", "w")); os.replace(out + ".tmp", out)
        el = time.time() - t0
        print(f"progress {tot['shots']}/{shots} ({el/60:.1f} min): accept exact {tot['accept_exact']:.4f} proxy {tot['accept_proxy']:.4f}; "
              f"eps exact {tot['eps_exact']:.3e} proxy {tot['eps_proxy']:.3e} ({tot['proxy_fail']} proxy fails); "
              f"triv: exact {tot['eps_exact_triv']:.3e} proxy {tot['eps_proxy_triv']:.3e}", flush=True)
    return tot


if __name__ == "__main__":
    if sys.argv[1] == "selftest":
        selftest()
    else:
        d, p, shots, out = int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
        seed = int(sys.argv[6]) if len(sys.argv) > 6 else 0
        batch = int(sys.argv[7]) if len(sys.argv) > 7 else 20000
        cfg = sys.argv[8] if len(sys.argv) > 8 else "serial"
        run(d, p, shots, out, seed, batch, cfg=cfg)
