"""Qudits at a given output rate (2026-11-08): the immediate-decision recipe of record (MEAS_DISTILL §5k) with F factories.
Factory = cultivation then escape, serial, on the escape's 702 qudits (the graft circuit's qudit count; cultivation runs in its
colour region).  Time per kept input = A (E[T_cult] + a_cult E[T_esc]) with early abort at the first firing postselected
detector (runs/surface/attempt_time.jsonl, ASAP ticks; 10 ticks = one surface-code round).  One output = 15.0 kept inputs
(the footprint's attempts/output / A) and one [[7,1,3]] s=1 round of 31 rounds (host_census schedule, tau_F = 0).
Qudits = 8 logical x 2 D^2 + F x 702; routing between factories and the round is not counted.
python scratchpad/rate_footprint.py"""
import json
import math

TICKS_PER_ROUND = 10
GADGET_ROUNDS = 31
FACTORY_Q = 702
# p: (cut, A = cultivation attempts per kept escape, attempts per output, D for 1e-10)
RECIPE = {1e-3: (7, 34.4, 517.2, 25), 7e-4: (7, 11.4, 171.2, 19), 5e-4: (7, 5.6, 84.0, 17)}
T = {r["p"]: r for r in map(json.loads, open("runs/surface/attempt_time.jsonl"))}
PERIODS = (31, 100, 300, 1000, 3000, 10000)
for p, (cut, A, att_out, D) in RECIPE.items():
    t = T[p]
    t_in = A * (t["cult_mean_ticks"] + t["cult_accept"] * t["esc_mean_ticks"]) / TICKS_PER_ROUND
    per_out = att_out / A * t_in                       # factory-rounds per output
    base = 8 * 2 * D * D
    print(f"p={p:.0e} cut {cut} D={D}: {t_in:.0f} factory-rounds per kept input, {per_out:.0f} per output; "
          f"base {base} q; one factory -> one output per {max(per_out, GADGET_ROUNDS):.0f} rounds ({base + FACTORY_Q} q)")
    row = []
    for P in PERIODS:
        F = math.ceil(per_out / P)
        row.append(f"{P}: F={F} {base + F * FACTORY_Q} q")
    print("   period (rounds/output) -> " + " | ".join(row))
    room = (10_000 - base) // FACTORY_Q
    print(f"   under 10k qudits: F={room}" + (f" -> one output per {max(per_out / room, GADGET_ROUNDS):.0f} rounds" if room > 0 else " (none)"))
