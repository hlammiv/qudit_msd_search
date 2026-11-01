"""Host-code census of the [[7,1,3]] s = 1 check gadget (MEAS_DISTILL §5j).

The planner charges the Clifford layer a_C p_C with a_C = 0.357 from the depolarizing census (cliff_census.py: every Pauli
of every gate equally likely) and p_C = one host logical error.  On the surface-code host the logical faults are not
depolarizing: a logical CSUM mostly fails by a single-wire logical Pauli, idling patches fail per round, a Fourier gate costs
a patch rotation, and a readout fails like one round.  This module schedules the gadget's circuit (meas_circuit.build) into
rounds (CSUM = 1 round, F/FD = tau_F rounds, R/M = 1 round), injects every single-wire logical Pauli Z^a / X^a (a = 1..d-1)
at every CSUM wire, idle round, Fourier gate, reset and readout, and counts the undetected bad outputs per unit of
per-location error ("exposure" N_cat, outputs per run per unit probability).  With the measured per-round and per-CSUM
logical error rates of the Z_4 host (footprint.BASE2_Z4 / CSUM_Z4 X-basis fits -> Z-type logical errors, the only type that
survives the check), the Clifford floor per output is

    floor(p, D) = N_gate q_gate(D) + (N_idle + tau_F N_F + N_meas + N_prep) q_round(D)

which replaces a_C p_L(D).  `python -m primepower_msd.host_census [d] [tau_F]` prints the census and the floor table; the
depolarizing census of cliff_census.py (0.357 per output at d = 4) is reproduced as a check.
"""
from __future__ import annotations

import itertools as it
import math
import sys
from collections import Counter, defaultdict

import numpy as np

from .meas_circuit import build
from .meas_distill import signed_rm
from .meas_sim import Protocol
from .zdsim import FrameSim

# X-basis (Z-type logical error) fits, d = 4, runs/surface/csum_d4.jsonl and zd_d4.jsonl (MEAS_DISTILL §5i, §5g)
BASE2_X = dict(A=0.193, p_th=6.95e-3)      # two patches, 2D rounds, no gate
CSUM_X = dict(A=0.214, p_th=6.90e-3)       # two patches, 2D rounds, one transversal CSUM, Pauli-frame decoder
HOST_X = dict(A=0.026, p_th=7.6e-3)        # single patch, D rounds, X basis (half of the both-basis fit A 0.053)


def _scale(fit, p, D):
    return fit["A"] * (p / fit["p_th"]) ** ((D + 1) / 2)


def q_round(p, D):
    """Z-type logical error per patch per round."""
    return _scale(BASE2_X, p, D) / (4 * D)


def q_gate(p, D):
    """Z-type logical error per wire per logical CSUM (its excess over idling, split over the two wires, plus the round)."""
    return max(_scale(CSUM_X, p, D) - _scale(BASE2_X, p, D), 0.0) / 2 + q_round(p, D)


def gadget(d):
    S, tau = signed_rm(1, 3)
    P = Protocol(S, tau, d, np.array([[0]]))
    C, inputs = build(P, 0.0, 0.0, 1e-3)
    return P, C, inputs


def schedule(C, tau_F=0):
    """ASAP rounds: returns (start time per op, free time per qudit per op, T_end).  Noise ops take no time."""
    DUR = {"CSUM": 1, "F": tau_F, "FD": tau_F, "R": 1, "M": 1}
    free = defaultdict(int)
    start = []
    for op in C.ops:
        name = op[0]
        if name == "CSUM":
            qs = [q for pair in op[1] for q in pair]
        elif name in ("F", "FD", "R", "M"):
            qs = list(op[1])
        else:
            start.append(None)
            continue
        s = max(free[q] for q in qs)
        start.append(s)
        for q in qs:
            free[q] = s + DUR[name]
    return start, dict(free), max(free.values())


def faults_logical(C, inputs, d, tau_F=0):
    """[(op index k, qudit, (x, z), weight, category)]: single-wire logical Paulis injected after op k."""
    start, _, T_end = schedule(C, tau_F)
    units = range(1, d)
    F = []
    last = {}                                      # qudit -> (op index of its last gate, free time after it)
    for k, op in enumerate(C.ops):
        name = op[0]
        if name == "CSUM":
            for (c, t) in op[1]:
                for q in (c, t):
                    for a in units:
                        F.append((k, q, (0, a), 1 / (d - 1), "gate")); F.append((k, q, (a, 0), 1 / (d - 1), "gate_x"))
                    _idle(F, last, q, start[k], k, d, units)
                    last[q] = (k, start[k] + 1)
        elif name in ("F", "FD"):
            for q in op[1]:
                _idle(F, last, q, start[k], k, d, units)
                for a in units:                    # the rotation: tau_F rounds of exposure on this wire
                    F.append((k, q, (0, a), tau_F / (d - 1), "fourier")); F.append((k, q, (a, 0), tau_F / (d - 1), "fourier_x"))
                last[q] = (k, start[k] + tau_F)
        elif name == "R":
            for q in op[1]:
                for a in units:                    # preparation error of a fresh patch: one round's worth
                    F.append((k, q, (a, 0), 1 / (d - 1), "prep")); F.append((k, q, (0, a), 1 / (d - 1), "prep_z"))
                last[q] = (k, start[k] + 1)
        elif name == "M":
            for q in op[1]:
                _idle(F, last, q, start[k], k, d, units)
                kk, _ = last[q]
                for a in units:                    # readout error: an X^a flip just before the measurement
                    F.append((kk, q, (a, 0), 1 / (d - 1), "meas"))
                last.pop(q, None)
    for q in list(last):                           # wires alive at the end (the outputs) idle until the gadget ends
        _idle(F, last, q, T_end, None, d, units)
    return F, T_end


def _idle(F, last, q, t_next, k_next, d, units):
    if q not in last:
        return
    k_prev, t_free = last[q]
    n = t_next - t_free
    if n > 0:
        for a in units:
            F.append((k_prev, q, (0, a), n / (d - 1), "idle")); F.append((k_prev, q, (a, 0), n / (d - 1), "idle_x"))


def faults_depolarizing(C, d):
    """cliff_census.py's fault set: every Pauli at every DEP1/DEP2/XERR site, equal weights (validation)."""
    F = []
    for k, op in enumerate(C.ops):
        name = op[0]
        if name in ("DEP1", "XERR"):
            vals = [(a, b) for a in range(d) for b in range(d) if (a, b) != (0, 0)] if name == "DEP1" else [(a, 0) for a in range(1, d)]
            for q in op[1]:
                for v in vals:
                    F.append((k, (q,), (v,), 1 / len(vals), name))
        elif name == "DEP2":
            vals = [v for v in it.product(range(d), repeat=4) if any(v)]
            for (c, t) in op[1]:
                for v in vals:
                    F.append((k, (c, t), ((v[0], v[1]), (v[2], v[3])), 1 / len(vals), name))
    return F


def census(C, inputs, d, faults, chunk=2000):
    """Per category: exposure N (sum of weight x bad outputs over undetected faults), Z-type only and any."""
    n = len(faults)
    det = np.zeros(n, bool); nbad = np.zeros(n, int); nbadZ = np.zeros(n, int)
    for lo in range(0, n, chunk):
        hi = min(n, lo + chunk); forced = {}
        for s in range(lo, hi):
            k, qs, pv = faults[s][0], faults[s][1], faults[s][2]
            qs = list(qs) if isinstance(qs, tuple) else [qs]
            pv = list(pv) if isinstance(pv[0], tuple) else [pv]
            mask = np.zeros(hi - lo, np.int64); mask[s - lo] = 1
            arr = np.zeros((hi - lo, len(qs), 2), np.int64)
            for i, (a, b) in enumerate(pv):
                arr[s - lo, i] = (a, b)
            forced.setdefault(k, []).append((mask, qs, arr))
        sim = FrameSim(C, hi - lo, np.random.default_rng(0), gauge=False).run(scale=0.0, forced=forced)
        det[lo:hi] = np.any(sim.detector_values() != 0, axis=1)
        z = sim.z[:, inputs] % d; x = sim.x[:, inputs] % d
        nbadZ[lo:hi] = np.count_nonzero(z, axis=1); nbad[lo:hi] = np.count_nonzero(z | x, axis=1)
    N = Counter(); NZ = Counter(); W = Counter()
    for s, f in enumerate(faults):
        w, cat = f[3], f[4]
        W[cat] += w
        if not det[s]:
            N[cat] += w * nbad[s]; NZ[cat] += w * nbadZ[s]
    return N, NZ, W


# Z-type exposures per output of the [[7,1,3]] s = 1 gadget at d = 4 (python -m primepower_msd.host_census 4 0): 45 CSUMs,
# depth 31 rounds; gate = CSUM-wire exposures, rounds = idle patch-rounds (Fourier gates, resets and readouts contribute no
# Z-type exposure; a tau_F = 15-round patch rotation per Fourier gate adds 14 idle rounds).
EXPOSURE_713 = dict(gate=4.667, rounds=5.333, rounds_tauF15=19.333)


# Fixed-D measurement of the excess (runs/surface/csum_excess.log, 1e6 shots): 2.0 +- 0.3 patch-rounds at D = 5, 3.4 +- 0.4 at
# D = 7, i.e. ~0.45 D -- the fitted 0.11 x 4D; a constant one extra round per patch is excluded (MEAS_DISTILL 5j).
def floor(p, D, NZ=None, tau_F=0, gate_model="fit"):
    """Clifford floor per output on the measured Z_4 host (Z-type exposures).  gate_model 'fit': q_gate from the CSUM-minus-
    memory fits (excess ~0.11 of the two-patch memory, i.e. ~0.4 D patch-rounds per CSUM); 'round': the CSUM layer costs one
    extra round per wire (q_gate = 2 q_round), the optimistic reading of the same data."""
    if NZ is None:
        NZ = dict(gate=EXPOSURE_713["gate"], idle=EXPOSURE_713["rounds_tauF15" if tau_F else "rounds"], fourier=0, meas=0, prep=0)
    rounds = NZ["idle"] + NZ["fourier"] + NZ["meas"] + NZ["prep"]      # fourier already carries tau_F in its weights
    qg = q_gate(p, D) if gate_model == "fit" else 2 * q_round(p, D)
    return NZ["gate"] * qg + rounds * q_round(p, D)


def main():
    d = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    tau_F = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    P, C, inputs = gadget(d)
    N, NZ, W = census(C, inputs, d, faults_depolarizing(C, d))
    print(f"depolarizing census (cliff_census check), d = {d}: per output any {sum(N.values()) / P.n_in:.3f}, Z-type {sum(NZ.values()) / P.n_in:.3f}")
    F, T_end = faults_logical(C, inputs, d, tau_F)
    start, free, _ = schedule(C, tau_F)
    n_csum = sum(1 for op in C.ops if op[0] == "CSUM")
    N, NZ, W = census(C, inputs, d, F)
    print(f"logical census, tau_F = {tau_F}: {n_csum} CSUMs, depth {T_end} rounds, {len(inputs)} output(s)")
    print("  category        exposure (weight)  undetected bad outputs: any   Z-type")
    for cat in ("gate", "gate_x", "idle", "idle_x", "fourier", "fourier_x", "meas", "prep", "prep_z"):
        if W[cat]:
            print(f"  {cat:14s} {W[cat]:8.2f}          {N[cat]:8.3f}  {NZ[cat]:8.3f}")
    print(f"  Z-type exposures per output: gate {NZ['gate'] / P.n_in:.3f}, rounds {(NZ['idle'] + NZ['fourier'] + NZ['meas'] + NZ['prep']) / P.n_in:.3f}")
    print("floor per output at p = 1e-3 (a_C p_L: 0.357 x 0.053 (p/7.55e-3)^((D+1)/2); csum: 0.357 x 0.216 (p/7.1e-3)^..):")
    for D in range(15, 29, 2):
        f = floor(1e-3, D, NZ, tau_F) / P.n_in; fr = floor(1e-3, D, NZ, tau_F, "round") / P.n_in
        old = 0.357 * 0.053 * (1e-3 / 7.55e-3) ** ((D + 1) / 2)
        cs = 0.357 * 0.216 * (1e-3 / 7.1e-3) ** ((D + 1) / 2)
        print(f"  D={D}: census {f:.2e} (gate = 1 extra round: {fr:.2e})  a_C p_L {old:.2e}  csum {cs:.2e}  (q_round {q_round(1e-3, D):.1e}, q_gate {q_gate(1e-3, D):.1e})")


if __name__ == "__main__":
    main()
