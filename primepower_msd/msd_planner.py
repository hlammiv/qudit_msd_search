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
#  Clifford fault locations per run L_C (meas_circuit CSUM encoder/decoder, slightly conservative count),
#  first-order Clifford coefficient a_C: per-output error a_C p_C from single Clifford faults (single-fault census of
#  the circuit, Z-type; equals the exact statevector value 0.359 for [[7,1,3]] s=1 at d = 4; A_CLIFF where no census))
STAGES_D4 = [
    ("O3 [[7,1,3]] s=1",          15.0,   15,     3, 3.148,  1,    7,   71,      0.36),
    ("O3 [[7,1,3]]",              29.0,   29,     3, 2.370,  1,    7,   140,     0.36),
    ("O3 [[15,7,3]]+K8",          9.571,  268,    3, 0.836,  28,   15,  2192,    1.02),
    ("O3 [[31,21,3]]+K22",        6.905,  1595,   3, 1.242,  231,  31,  26796, None),
    ("O3 [[63,51,3]]+K52",        5.941,  7878,   3, 2.465,  1326, 63,  273520, None),
    ("O3 [[127,113,3]]+K114",     5.496,  35397,  3, 5.076,  6441, 127, 2550408, None),
    ("O7 [[63,21,7]]+1470",       19.0,   27930,  7, 343.3,  1470, 63,  534450, None),
    ("O7 [[64,26,7]]+2080",       15.769, 32800,  7, 360.0,  2080, 64,  1013120, None),
    ("O7 [[112,50,7]]+9500",      14.44,  137180, 7, 1927.2, 9500, 112, 6314080, None),
]
# L_C = 2 N + M (2 net + 4 n - 3 k + 6) with net = CSUMs per dense encoder (sparse_encoder.py cuts net to
# 1384 / 4721 for [[64,26]] / [[112,50]], i.e. L_C 712640 / 5574220; order-3 stages 1.5-2.6x fewer, see MEAS_DISTILL 5e).
# Other dimensions: order-3 coefficients from runs/meas_distill/order3_d{d}.json (same T/out, L, L_C: the circuits
# are d-independent), order-7 coefficients from the exact counts below, a_C from the d = 8 census where available.
ORDER7_COEF = {  # (d, stage name) -> per-output coefficient c (runs/meas_distill/count_*_d{d}.log)
    (8, "O7 [[63,21,7]]+1470"): 2.126,
    (16, "O7 [[63,21,7]]+1470"): 0.0220,
    (8, "O7 [[64,26,7]]+2080"): 2.229,
    (16, "O7 [[64,26,7]]+2080"): 0.0230,
}
A_C_D8 = {"O3 [[7,1,3]] s=1": 0.147, "O3 [[7,1,3]]": 0.147, "O3 [[15,7,3]]+K8": 0.42}


def stages_for_d(d: int, path: str | None = None):
    """STAGES for d = 8, 16 (STAGES_D4 for d = 4): order-3 rows re-coefficiented from order3_d{d}.json, order-7 rows
    only where an exact count exists."""
    import json
    import os
    if d == 4:
        return STAGES_D4
    path = path or os.path.join(os.path.dirname(__file__), "..", "runs", "meas_distill", f"order3_d{d}.json")
    o3 = json.load(open(path))
    out = []
    for name, t, L, W, c, inputs, n, LC, aC in STAGES_D4:
        if W == 3:
            m = {7: "3", 15: "4", 31: "5", 63: "6", 127: "7"}[n]
            if m not in o3:
                continue
            if "s=1" in name:
                c_new = {8: 0.539, 16: 0.1141}.get(d)        # single-check [[7,1,3]]: separate exact count
                if c_new is None:
                    continue
            else:
                c_new = o3[m]["coef"]
            out.append((name, t, L, W, c_new, inputs, n, LC, (A_C_D8 if d == 8 else {}).get(name)))
        elif (d, name) in ORDER7_COEF:
            out.append((name, t, L, W, ORDER7_COEF[(d, name)], inputs, n, LC, None))
    return out


# First-order Clifford term: per-output error a * p_C.  The Pauli-frame circuit model gives a ~ 3-5 when X-type
# output errors are counted and 0.36-1.0 for Z-type; the exact statevector simulation of the [[7,1,3]] gadget
# (real T gates, d = 4) gives 0.359 = the Z-type value: X-type faults are projected out or detected by the check.
# A_CLIFF is the conservative fallback for stages without a census (a_C None).
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
            for name, t, L, W, c, _, _, LC, aC in seq:
                if L * p + LC * pC > valid_Lp:
                    ok = False
                    break
                q = c * p ** W + (aC if aC is not None else a_cliff) * pC
                if q >= p:
                    ok = False
                    break
                cost *= t * math.exp(L * p + LC * pC)
                trace.append((name, p, q, math.exp(-L * p - LC * pC)))
                p = q
            if ok and p <= target and (best is None or cost < best[0]):
                best = (cost, trace)
    return best


def table(ps=(1e-3, 3e-4, 1e-4, 1e-5, 1e-6), targets=(1e-10, 1e-15, 1e-20), stages=STAGES_D4):
    rows = []
    for p0 in ps:
        for tg in targets:
            b = plan(p0, tg, stages=stages)
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
    "grow3->5 (56 q)": {1e-3: (0.452, 9.07e-7), 5e-4: (0.672, 1.49e-7)},  # 1e8 shots, sweep_grow35_d4_long
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
    import sys
    d = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    print(f"d = {d}")
    for p0, tg, cost, trace in table(stages=stages_for_d(d)):
        if cost is None:
            print(f"p={p0:.0e} target={tg:.0e}: no plan")
            continue
        steps = " -> ".join(f"{n} (p_in {p:.1e}, out {q:.1e}, acc {a:.2f})" for n, p, q, a in trace)
        print(f"p={p0:.0e} target={tg:.0e}: {cost:8.1f} T/output   {steps}")
