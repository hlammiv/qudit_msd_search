"""Vectorised Pauli-frame simulator for Clifford circuits on Z_d qudits (any d), stim-style.

Frame = X^x Z^z per qudit (phases ignored).  Conjugation rules (verified against dense matrices in
tests/test_zdsim.py):
    F   (F|j> = sum_k w^{jk}|k>/sqrt d):  (x, z) -> (-z, x)
    FD  (= F^dag):                         (x, z) -> (z, -x)
    S^a (diag zeta_{2d}^{a x^2}, even d):  z += a x
    CSUM(c -> t, w)  |a,b> -> |a, b + w a>: x_t += w x_c,  z_c -= w z_t
    M   (Z-basis measurement):             record x (outcome shift); then z += uniform (gauge)
    R   (reset to |0>):                    x = z = 0; z = uniform (gauge)
    FF  (classically controlled Pauli):    x_q / z_q += coef * m_i  (feed-forward on a recorded outcome; the
                                           reference run has all outcomes 0, so only the flip enters the frame)
Gauge randomisation (stim's trick) makes every non-deterministic measurement uniformly random, so a declared
detector that fires in a noiseless run is NOT deterministic -> `check_determinism`.

Noise: DEP1(p) uniform non-identity single-qudit Pauli; DEP2(p) uniform non-identity two-qudit Pauli;
XERR(p) random X^j (j != 0) — used after resets and before measurements (= measurement error).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Circuit:
    d: int
    nq: int = 0
    ops: list = field(default_factory=list)
    n_meas: int = 0
    detectors: list = field(default_factory=list)       # list of [(meas_index, coeff), ...]
    det_tags: list = field(default_factory=list)

    def q(self, k: int = 1) -> list[int]:
        out = list(range(self.nq, self.nq + k)); self.nq += k; return out

    def _add(self, *op):
        self.ops.append(op)

    # gates
    def R(self, qs): self._add("R", list(qs))
    def F(self, qs): self._add("F", list(qs))
    def FD(self, qs): self._add("FD", list(qs))
    def S(self, qs, a=1):
        qs = list(qs); a = np.broadcast_to(np.asarray(a, dtype=np.int64), (len(qs),))
        self._add("S", qs, [int(v) for v in a])
    def CSUM(self, pairs, w=1):
        pairs = list(pairs); w = np.broadcast_to(np.asarray(w, dtype=np.int64), (len(pairs),))
        self._add("CSUM", pairs, [int(v) for v in w])
    def M(self, qs) -> list[int]:
        qs = list(qs); idx = list(range(self.n_meas, self.n_meas + len(qs))); self.n_meas += len(qs)
        self._add("M", qs, idx); return idx
    def FF(self, terms):
        """terms: [(meas_index, qudit, 'x'|'z', coef)] -- Pauli feed-forward (gauge fixing / frame updates)."""
        self._add("FF", [(int(i), int(q), ax, int(c)) for i, q, ax, c in terms])
    # noise
    def DEP1(self, qs, p): self._add("DEP1", list(qs), p)
    def DEP2(self, pairs, p): self._add("DEP2", list(pairs), p)
    def XERR(self, qs, p): self._add("XERR", list(qs), p)
    def ZERR(self, qs, p): self._add("ZERR", list(qs), p)

    def detector(self, terms, tag=""):
        self.detectors.append([(int(i), int(c)) for i, c in terms]); self.det_tags.append(tag)

    def noise_locations(self) -> int:
        return sum(1 for op in self.ops if op[0] in ("DEP1", "DEP2", "XERR", "ZERR"))


class FrameSim:
    def __init__(self, circ: Circuit, shots: int, rng=None, gauge: bool = True):
        self.c, self.shots, self.d = circ, shots, circ.d
        self.rng = rng if rng is not None else np.random.default_rng()
        self.gauge = gauge
        self.x = np.zeros((shots, circ.nq), dtype=np.int64)
        self.z = np.zeros((shots, circ.nq), dtype=np.int64)
        self.meas = np.zeros((shots, circ.n_meas), dtype=np.int64)

    def _rand_nonid(self, k, size):
        """uniform non-zero vectors in Z_d^k, shape (size, k)."""
        r = self.rng.integers(1, self.d ** k, size=size)
        return np.stack([(r // self.d ** j) % self.d for j in range(k)], axis=1)

    def run(self, scale: float = 1.0, forced=None):
        """forced: optional dict op_index -> [(shot_mask, qudits, pauli_array)], faults applied after op k."""
        d, x, z = self.d, self.x, self.z
        for k, op in enumerate(self.c.ops):
            name = op[0]
            if name == "R":
                q = op[1]; x[:, q] = 0
                z[:, q] = self.rng.integers(0, d, size=(self.shots, len(q))) if self.gauge else 0
            elif name == "F":
                q = op[1]; x[:, q], z[:, q] = (-z[:, q]) % d, x[:, q].copy()
            elif name == "FD":
                q = op[1]; x[:, q], z[:, q] = z[:, q].copy(), (-x[:, q]) % d
            elif name == "S":
                for qq, a in zip(op[1], op[2]):
                    z[:, qq] = (z[:, qq] + a * x[:, qq]) % d
            elif name == "CSUM":
                for (cq, tq), w in zip(op[1], op[2]):
                    x[:, tq] = (x[:, tq] + w * x[:, cq]) % d
                    z[:, cq] = (z[:, cq] - w * z[:, tq]) % d
            elif name == "M":
                q, idx = op[1], op[2]
                self.meas[:, idx] = x[:, q]
                if self.gauge:
                    z[:, q] = (z[:, q] + self.rng.integers(0, d, size=(self.shots, len(q)))) % d
            elif name == "FF":
                for i, qq, ax, cf in op[1]:
                    arr = x if ax == "x" else z
                    arr[:, qq] = (arr[:, qq] + cf * self.meas[:, i]) % d
            elif name in ("DEP1", "XERR", "ZERR"):
                q, p = op[1], op[2] * scale
                if p > 0:
                    hit = self.rng.random((self.shots, len(q))) < p
                    if name == "DEP1":
                        P = self._rand_nonid(2, hit.shape).reshape(*hit.shape, 2)
                        x[:, q] = (x[:, q] + hit * P[..., 0]) % d; z[:, q] = (z[:, q] + hit * P[..., 1]) % d
                    else:
                        v = self.rng.integers(1, d, size=hit.shape)
                        if name == "XERR": x[:, q] = (x[:, q] + hit * v) % d
                        else: z[:, q] = (z[:, q] + hit * v) % d
            elif name == "DEP2":
                pairs, p = op[1], op[2] * scale
                if p > 0:
                    hit = self.rng.random((self.shots, len(pairs))) < p
                    P = self._rand_nonid(4, hit.shape).reshape(*hit.shape, 4)
                    a = [pp[0] for pp in pairs]; b = [pp[1] for pp in pairs]
                    x[:, a] = (x[:, a] + hit * P[..., 0]) % d; z[:, a] = (z[:, a] + hit * P[..., 1]) % d
                    x[:, b] = (x[:, b] + hit * P[..., 2]) % d; z[:, b] = (z[:, b] + hit * P[..., 3]) % d
            if forced is not None and k in forced:
                for mask, qs, P in forced[k]:      # P: (shots, len(qs), 2)
                    x[:, qs] = (x[:, qs] + mask[:, None] * P[..., 0]) % d
                    z[:, qs] = (z[:, qs] + mask[:, None] * P[..., 1]) % d
        return self

    def detector_values(self) -> np.ndarray:
        D = np.zeros((self.shots, len(self.c.detectors)), dtype=np.int64)
        for j, terms in enumerate(self.c.detectors):
            for i, cf in terms:
                D[:, j] += cf * self.meas[:, i]
        return D % self.d


def add_tick_idle(circ: Circuit, p_idle: float) -> Circuit:
    """Copy of `circ` with idle noise at every tick: the gates are scheduled as soon as possible in circuit order (one tick
    per gate; the pairs of one CSUM op take one tick each, in order, so a fan-out serialises on its shared qudit); a qudit
    is live from its reset until it is measured, and every tick a live qudit waits gets DEP1(p_idle) -- k waiting ticks
    compose exactly into one DEP1 inserted before the gate it waits for; live qudits also wait until the circuit's last
    tick.  Resets are just in time: a freshly reset qudit (and its single-qudit preparation gates) waits for nothing before
    its first two-qudit gate.  Noise, detector and feed-forward ops take no time.  Measurement indices and detectors are
    unchanged."""
    if p_idle <= 0:
        return circ
    d = circ.d
    lam = p_idle * d * d / (d * d - 1)                    # DEP1(q) = (1 - lam) rho + lam (fully depolarised)

    def comp(k):
        return (1 - (1 - lam) ** k) * (d * d - 1) / (d * d)

    out = Circuit(d, circ.nq, [], circ.n_meas, circ.detectors, circ.det_tags)
    t = [0] * circ.nq
    live = [False] * circ.nq
    fresh = [False] * circ.nq                             # reset, not yet in a two-qudit gate: prepared just in time

    def emit_idle(waits):
        by_k = {}
        for q, k in waits.items():
            if k > 0:
                by_k.setdefault(k, []).append(q)
        for k in sorted(by_k):
            out.ops.append(("DEP1", by_k[k], comp(k)))

    for op in circ.ops:
        name = op[0]
        if name in ("F", "FD", "S", "R", "M"):
            for q in op[1]:
                t[q] += 1                                # single-qudit gates run as soon as the qudit is free
            out.ops.append(op)
            if name == "R":
                for q in op[1]:
                    live[q] = fresh[q] = True
            elif name == "M":
                for q in op[1]:
                    live[q] = fresh[q] = False
        elif name == "CSUM":
            waits = {}
            for (a, b) in op[1]:
                s0 = max(t[a], t[b])
                for q in (a, b):
                    if live[q] and not fresh[q] and s0 > t[q]:
                        waits[q] = waits.get(q, 0) + s0 - t[q]
                    fresh[q] = False
                t[a] = t[b] = s0 + 1
            emit_idle(waits)
            out.ops.append(op)
        else:
            out.ops.append(op)
    T = max(t) if t else 0
    emit_idle({q: T - t[q] for q in range(circ.nq) if live[q]})
    return out


def check_determinism(circ: Circuit, shots: int = 256, seed: int = 0) -> list[int]:
    """Indices of declared detectors that fire in noiseless gauge-randomised runs (should be [])."""
    s = FrameSim(circ, shots, np.random.default_rng(seed), gauge=True).run(scale=0.0)
    return [int(j) for j in np.nonzero(np.any(s.detector_values(), axis=0))[0]]
