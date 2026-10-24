"""Physical-qudit footprint of cultivation -> HHPW distillation on a topological host code (d = 2^k ring qudits).

Model (every number here is an assumption to revisit, not a measurement):
* Host code: logical error per logical gate p_L(D) = A (p / p_th)^((D+1)/2), footprint F(D) physical qudits per
  logical qudit (surface-code-like: 2 D^2).  Defaults A = 0.1, p_th = 5e-3 (qubit surface-code numbers; the ring
  Z_{2^k} surface code has no measured threshold yet).
* Cultivation (msd_planner.CULT_D4) delivers states at p_cult with acceptance a_cult on a 56-qudit block; moving one
  into the host code ("escape") costs one logical-gate error p_L(D).
* An HHPW stage on logical qudits: per-output error c p^W + a_C p_C with p_C = p_L(D) (a_C per stage from the
  single-fault census, exact 0.36 for [[7,1,3]]; A_CLIFF = 5 where unknown), acceptance
  exp(-L p - L_C p_C) (msd_planner.plan); logical qudits in flight = inputs + (n - k) block + 1 control.
Output: for each (p, target) the smallest host distance that reaches the target and the cheapest pipeline at that
distance, with physical qudits = logical qudits x F(D) + 56 (one cultivation factory) and cultivation attempts per
delivered output.
"""
from __future__ import annotations

import itertools as it
import math
import re
from dataclasses import dataclass

from .msd_planner import A_CLIFF, CULT_D4, STAGES_D4, plan

CULT_KEY = "grow3->5 (56 q)"

# Measured on the signed Z_4 rotated surface code with the 2-adic matching decoder (surface_zd, runs/surface/zd_d4.jsonl):
# circuit-level threshold ~0.8% (D=7/9 and 9/11 crossings 7.7-8.3e-3, same at d=8), fit of p_fail(D rounds) =
# A (p/p_th)^((D+1)/2) over the 129 sub-threshold points of d=4 and d=8, both bases: A ~ 0.053, p_th ~ 0.0076.
HOST_Z4 = dict(p_th=7.55e-3, A=0.053, name="Z4/Z8 rotated surface (measured)")
CULT_QUDITS = 56


@dataclass(frozen=True)
class Host:
    p_th: float = 5e-3
    A: float = 0.1
    name: str = "surface-like"

    def p_logical(self, p: float, D: int) -> float:
        return self.A * (p / self.p_th) ** ((D + 1) / 2)

    def qudits(self, D: int) -> int:
        return 2 * D * D


def logical_qudits(stage) -> int:
    """inputs + (n - k) reusable block qudits + 1 control (k read from the stage name [[n,k,D]])."""
    name, t, L, W, c, inputs, n, LC, aC = stage
    k = int(re.search(r"\[\[\d+,(\d+)", name).group(1))
    return inputs + (n - k) + 1


def pipelines(p: float, target: float, host: Host = Host(), stages=STAGES_D4, cult=CULT_D4, a_cliff: float = A_CLIFF,
              max_levels: int = 2, Dmax: int = 41, valid_Lp: float = 3.0):
    """Every (host distance, stage sequence) reaching `target` at physical error p: dicts with D, pipeline, error,
    logical and physical qudits, T per accepted output, cultivation attempts per output."""
    a_cult, p_cult = cult[CULT_KEY][p]
    out = []
    for D in range(3, Dmax + 1, 2):
        pL = host.p_logical(p, D)
        p_in = p_cult + pL                                   # escape into the host code
        if p_in <= target:                                   # cultivation only
            out.append(dict(D=D, pipeline="cultivation only", error=p_in, logical=1, physical=host.qudits(D) + CULT_QUDITS,
                            cult_attempts=1 / a_cult, T_per_out=1.0))
        if min((s[8] if s[8] is not None else a_cliff) for s in stages) * pL >= target:
            continue                                         # first-order Clifford floor: no stage reaches the target
        for depth in range(1, max_levels + 1):
            for seq in it.product(stages, repeat=depth):
                q, cost, ok = p_in, 1.0, True
                for name, t, L, W, c, _, _, LC, aC in seq:
                    if L * q + LC * pL > valid_Lp:
                        ok = False
                        break
                    q_out = c * q ** W + (aC if aC is not None else a_cliff) * pL
                    if q_out >= q:
                        ok = False
                        break
                    cost *= t * math.exp(L * q + LC * pL)
                    q = q_out
                if ok and q <= target:
                    lq = max(logical_qudits(s) for s in seq)
                    out.append(dict(D=D, pipeline=" -> ".join(s[0] for s in seq), error=q, logical=lq,
                                    physical=lq * host.qudits(D) + CULT_QUDITS, cult_attempts=cost / a_cult, T_per_out=cost))
    return out


def cheapest(p, target, key="physical", **kw):
    cands = pipelines(p, target, **kw)
    return min(cands, key=lambda r: (r[key], r["cult_attempts"])) if cands else None


def table(ps=(1e-3, 5e-4), targets=(1e-10, 1e-12, 1e-15), **kw):
    rows = []
    for p in ps:
        for tg in targets:
            rows.append((p, tg, cheapest(p, tg, **kw), cheapest(p, tg, key="cult_attempts", **kw)))
    return rows


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "measured":
        host = Host(**HOST_Z4)
    else:
        host = Host(p_th=float(sys.argv[1]), A=float(sys.argv[2]) if len(sys.argv) > 2 else 0.1) if len(sys.argv) > 1 else Host()
    print(f"host: p_th {host.p_th:.0e}, A {host.A}, F(D) = 2 D^2; cultivation {CULT_KEY}; a_cliff {A_CLIFF}")
    for p, tg, r, r2 in table(host=host):
        if r is None:
            print(f"p={p:.0e} target={tg:.0e}: unreachable (D <= 41)")
            continue
        print(f"p={p:.0e} target={tg:.0e}: fewest qudits D={r['D']} {r['pipeline']}: {r['physical']} physical "
              f"({r['logical']} logical), err {r['error']:.1e}, {r['cult_attempts']:.1f} cultivation attempts/output"
              + (f" | fewest attempts: D={r2['D']} {r2['pipeline']} {r2['physical']} q, {r2['cult_attempts']:.1f} att/out"
                 if r2 != r else ""))
