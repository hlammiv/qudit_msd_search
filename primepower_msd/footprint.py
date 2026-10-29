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

# Measured escape (graft_zd / graft_run variants, runs/surface/graft_variants.jsonl; MEAS_DISTILL §5h): colour-code d=5 state
# grafted into the D=15 surface code with Gidney's d1=5 -> d2=15 graft (229 data qudits), 5 grafted + 5 final rounds,
# colour-region plaquettes postselected during the grafted rounds (joint over both sectors) plus decoder-certified
# first-order postselection.  Per physical error rate p: (escape error per kept attempt = X-sector + Z-sector failure
# rates; joint discard fraction per attempt).  The error is a one-off cost that adds to p_cult; it does not scale with the
# host distance (the D=15 patch is fixed by the graft).  2e5 attempts per sector and point.
ESCAPE_Z4 = {1e-3: (1.8e-4, 0.43), 5e-4: (2.5e-5, 0.25)}       # one colour class per sector (X: 2, Z: 0); X 9.7e-5 + Z 7.9e-5 at 1e-3
ESCAPE_Z4_ALL = {1e-3: (6.6e-5, 0.62), 5e-4: (1.0e-5, 0.39)}   # all three classes, both sectors; X 4.0e-5 + Z 2.7e-5 (5 events) at 1e-3
# one class + decoder-certified second-order postselection on the colour region (13M/16M fault pairs decoded at build time),
# first-order cost of the certified patterns capped at 0.2% (degenerate_cap=2): the best recipe measured
ESCAPE_Z4_PAIRS = {1e-3: (9.1e-5, 0.43), 5e-4: (1.3e-5, 0.24)}  # X 6.7e-5 + Z 2.5e-5 (52 events in 2e6 attempts) at 1e-3; 2 events at 5e-4
# ... plus complementary-gap postselection (PeelDecoder(gap=True), cut at 13): no failure in 1.14M kept attempts at 1e-3 ->
# 90% CL upper bound 4e-6 per sector; the recipe of record (MEAS_DISTILL 5h).  At 5e-4: no failure in 0.75M kept attempts
# per sector (1e6 attempts, 11 events per sector before the cut) -> < 3e-6 per sector, 24.9% discard.
ESCAPE_Z4_GAP = {1e-3: (8e-6, 0.437), 5e-4: (6e-6, 0.249)}

# Measured logical CSUM (logical_csum.py, runs/surface/csum_d4.jsonl; MEAS_DISTILL §5i): transversal CSUM between two
# D-patches with D rounds before and after, failure of either logical, Pauli-frame (sequential) decoder; fit
# p_csum(D) = A (p/p_th)^((D+1)/2).  BASE2_Z4: the same two patches with no gate (two-patch memory over 2D rounds).
# Fits over 17 sub-threshold points per basis (D = 3..9, p = 1e-3..6e-3, 1e5/4e4 shots): CSUM Z A 0.219 p_th 7.34e-3,
# X A 0.214 p_th 6.90e-3; two-patch memory Z A 0.190 p_th 7.34e-3, X A 0.193 p_th 6.95e-3 (basis averages below).  The gate
# costs 1.1-1.2x the idle two-patch memory at the same threshold; the gate's excess over idling is ~0.5 p_L(D) of one patch.
CSUM_Z4 = dict(p_th=7.1e-3, A=0.216, name="Z4 logical CSUM, Pauli-frame decoder (measured)")
BASE2_Z4 = dict(p_th=7.15e-3, A=0.19, name="Z4 two-patch memory, 2D rounds (measured)")


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
              max_levels: int = 2, Dmax: int = 41, valid_Lp: float = 3.0, escape=None, clifford=None):
    """Every (host distance, stage sequence) reaching `target` at physical error p: dicts with D, pipeline, error,
    logical and physical qudits, T per accepted output, cultivation attempts per output.  escape: {p: (p_escape, discard)}
    measured escape (default: one host logical error p_L(D), no discard).  clifford: Host whose p_logical(p, D) is the
    error per logical CSUM (the a_C and L_C fault locations are CSUMs); default one host logical error p_L(D)."""
    a_cult, p_cult = cult[CULT_KEY][p]
    if escape is not None:
        p_esc, disc = escape[p]
        a_cult = a_cult * (1 - disc)                         # an escape attempt that is discarded costs a cultivation
    out = []
    for D in range(3, Dmax + 1, 2):
        pL = host.p_logical(p, D)
        pC = pL if clifford is None else clifford.p_logical(p, D)   # error per logical CSUM
        p_in = p_cult + (pL if escape is None else p_esc)   # escape into the host code
        if p_in <= target:                                   # cultivation only
            out.append(dict(D=D, pipeline="cultivation only", error=p_in, logical=1, physical=host.qudits(D) + CULT_QUDITS,
                            cult_attempts=1 / a_cult, T_per_out=1.0))
        if min((s[8] if s[8] is not None else a_cliff) for s in stages) * pC >= target:
            continue                                         # first-order Clifford floor: no stage reaches the target
        for depth in range(1, max_levels + 1):
            for seq in it.product(stages, repeat=depth):
                q, cost, ok = p_in, 1.0, True
                for name, t, L, W, c, _, _, LC, aC in seq:
                    if L * q + LC * pC > valid_Lp:
                        ok = False
                        break
                    q_out = c * q ** W + (aC if aC is not None else a_cliff) * pC
                    if q_out >= q:
                        ok = False
                        break
                    cost *= t * math.exp(L * q + LC * pC)
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
    escape = None
    clifford = None
    if "csum" in sys.argv[2:]:                       # conservative: the whole two-patch CSUM experiment error per CSUM location
        clifford = Host(**CSUM_Z4)
        sys.argv.remove("csum")
    if "csum-excess" in sys.argv[2:]:                # the gate's measured excess over idling, on top of one p_L(D) per location
        class _Excess(Host):
            def p_logical(self, p, D):
                return Host(**HOST_Z4).p_logical(p, D) + Host(**CSUM_Z4).p_logical(p, D) - Host(**BASE2_Z4).p_logical(p, D)
        clifford = _Excess(**{**CSUM_Z4, "name": "p_L(D) + measured CSUM excess"})
        sys.argv.remove("csum-excess")
    if len(sys.argv) > 1 and sys.argv[1] in ("measured", "escape", "escape-all", "escape-pairs", "escape-gap"):
        host = Host(**HOST_Z4)
        escape = {"escape": ESCAPE_Z4, "escape-all": ESCAPE_Z4_ALL, "escape-pairs": ESCAPE_Z4_PAIRS,
                  "escape-gap": ESCAPE_Z4_GAP}.get(sys.argv[1])
    else:
        host = Host(p_th=float(sys.argv[1]), A=float(sys.argv[2]) if len(sys.argv) > 2 else 0.1) if len(sys.argv) > 1 else Host()
    print(f"host: p_th {host.p_th:.0e}, A {host.A}, F(D) = 2 D^2; cultivation {CULT_KEY}; a_cliff {A_CLIFF}; "
          f"escape {'measured ' + str(escape) if escape else 'one p_L(D)'}; "
          f"CSUM {clifford.name if clifford else 'one p_L(D)'}")
    ps = tuple(sorted(escape, reverse=True)) if escape else (1e-3, 5e-4)
    for p, tg, r, r2 in table(ps=ps, host=host, escape=escape, clifford=clifford):
        if r is None:
            print(f"p={p:.0e} target={tg:.0e}: unreachable (D <= 41)")
            continue
        print(f"p={p:.0e} target={tg:.0e}: fewest qudits D={r['D']} {r['pipeline']}: {r['physical']} physical "
              f"({r['logical']} logical), err {r['error']:.1e}, {r['cult_attempts']:.1f} cultivation attempts/output"
              + (f" | fewest attempts: D={r2['D']} {r2['pipeline']} {r2['physical']} q, {r2['cult_attempts']:.1f} att/out"
                 if r2 != r else ""))
