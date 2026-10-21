"""Concatenated HHPW distillation planner over Z_{2^k} (full post-selection per stage).

Each stage: T per output t, fault locations per run L (every single fault is detected => acceptance exp(-L p)),
per-output output error c * p^W (exact leading-order counts from meas_sim, d = 4, uniform nonzero fault values;
inputs and the stage's own T gates both at the incoming error p).  A plan is a sequence of stages; the cost per
final output is prod_i t_i * exp(L_i p_i).
"""
from __future__ import annotations

import itertools as it
import math

# (name, T per output, T-fault locations per run, W, per-output coefficient c, inputs per run, inner n,
#  Clifford fault locations per run L_C (meas_circuit CSUM encoder/decoder, slightly conservative count))
STAGES_D4 = [
    ("O3 [[7,1,3]]",              29.0,   29,     3, 2.370,  1,    7,   71),
    ("O3 [[15,7,3]]+K8",          9.571,  268,    3, 0.836,  28,   15,  2192),
    ("O3 [[31,21,3]]+K22",        6.905,  1595,   3, 1.242,  231,  31,  26796),
    ("O3 [[63,51,3]]+K52",        5.941,  7878,   3, 2.465,  1326, 63,  273520),
    ("O3 [[127,113,3]]+K114",     5.496,  35397,  3, 5.076,  6441, 127, 2550408),
    ("O7 [[63,21,7]]+1470",       19.0,   27930,  7, 343.3,  1470, 63,  534450),
    ("O7 [[64,26,7]]+2080",       15.769, 32800,  7, 360.0,  2080, 64,  1013120),
]
# First-order Clifford term: per-output error a * p_C (circuit-level, d = 4: a ~ 3-5 counting X and Z output errors,
# 0.4-1.1 for Z only; runs/meas_distill/circuit_clifford_d4.json).  A_CLIFF is the conservative value.
A_CLIFF = 5.0


def plan(p0: float, target: float, stages=STAGES_D4, max_levels: int = 3, valid_Lp: float = 3.0, pC: float = 0.0,
         a_cliff: float = A_CLIFF):
    """Cheapest stage sequence reaching per-output error <= target.  Stages with L p > valid_Lp are excluded
    (acceptance < 5 %, and the leading-order error formula is no longer trustworthy there).
    Clifford noise (circuit level): each stage adds a_cliff * pC to its per-output error (first order: faults in
    the final decoder land on outputs) and multiplies its acceptance by exp(-L_C pC)."""
    best = None
    for depth in range(1, max_levels + 1):
        for seq in it.product(stages, repeat=depth):
            p, cost, ok = p0, 1.0, True
            trace = []
            for name, t, L, W, c, _, _, LC in seq:
                if L * p + LC * pC > valid_Lp:
                    ok = False
                    break
                q = c * p ** W + a_cliff * pC
                if q >= p:
                    ok = False
                    break
                cost *= t * math.exp(L * p + LC * pC)
                trace.append((name, p, q, math.exp(-L * p - LC * pC)))
                p = q
            if ok and p <= target and (best is None or cost < best[0]):
                best = (cost, trace)
    return best


def table(ps=(1e-3, 3e-4, 1e-4, 1e-5, 1e-6), targets=(1e-10, 1e-15, 1e-20)):
    rows = []
    for p0 in ps:
        for tg in targets:
            b = plan(p0, tg)
            if b is None:
                rows.append((p0, tg, None, None))
                continue
            cost, trace = b
            rows.append((p0, tg, cost, trace))
    return rows


# Cultivation front end (d = 4, sweeps in runs/cultivation; p_L upper bound used when 0 failures were seen)
CULT_D4 = {
    "injection":    {1e-3: (1.0, 1e-3), 5e-4: (1.0, 5e-4)},
    "D3 (20 q)":    {1e-3: (0.727, 2.34e-5), 5e-4: (0.853, 5.75e-6)},
    "grow3->5 (56 q)": {1e-3: (0.452, 8.85e-7), 5e-4: (0.672, 4.5e-7)},
}


def end_to_end(p_phys: float, target: float, pC: float = 0.0):
    """For each front end: (front-end attempts per final output, distillation T per output, plan)."""
    out = {}
    for fe, tab in CULT_D4.items():
        acc, p0 = tab[p_phys]
        b = plan(p0, target, pC=pC) if p0 > target else (1.0, [])
        if b is None:
            out[fe] = None
            continue
        out[fe] = (b[0] / acc, b[0], b[1])
    return out


if __name__ == "__main__":
    for p0, tg, cost, trace in table():
        if cost is None:
            print(f"p={p0:.0e} target={tg:.0e}: no plan")
            continue
        steps = " -> ".join(f"{n} (p_in {p:.1e}, out {q:.1e}, acc {a:.2f})" for n, p, q, a in trace)
        print(f"p={p0:.0e} target={tg:.0e}: {cost:8.1f} T/output   {steps}")
