"""Proved protocol savings for the level-set codes (LEVEL_SETS.md) without any census: the weight-d
term is bounded by Singleton too (a d-set carries at most p-1 kernel vectors), so
L_d <= A_d <= C(n,d)(p-1), and every heavier accepted error is a failure (qmsd.rigorous_curve).
Baseline: the standard library at p plus the whole Reed-Solomon family (exact closed-form curves).
usage: levelset_cost.py OUT.json"""
import sys, json, math
from math import comb
from qmsd.objective import standard_library, DEFAULT_GRID
from qmsd.protocol import CodeCurve, best_protocol
from qmsd.rs_family import rs_profile, rs_parameters, r_max_rs
from qmsd.rigorous_curve import RigorousBoundCurve, rigorous_saving
codes = [("[[13332,1309,22]]_11", 11, 13332, 1309, 22), ("[[78455,5066,>=68]]_17", 17, 78455, 5066, 68), ("[[123284,7037,>=57]]_19", 19, 123284, 7037, 57)]
out = {}
for name, p, n, k, d in codes:
    base = standard_library(p)
    for kk in range(1, r_max_rs(p) + 1):
        nn, _, dd = rs_parameters(p, kk)
        if dd >= 2: base.append(CodeCurve(f"[[{nn},{kk},{dd}]]_{p} RS", p, nn, kk, rs_profile(p, kk)))
    A = comb(n, d) * (p - 1)
    c = RigorousBoundCurve(name, p, n, k, d, L_d=float(A), A_d=A)
    res = rigorous_saving(c, base)
    rows = []
    for r in res["rows"]:
        rows.append({kk: r[kk] for kk in ("eps0", "target", "base", "with_upper", "saving_proved", "uses_proved", "protocol_proved", "reaches_target_check")})
        b = best_protocol(base, r["eps0"], r["target"], max_depth=5)
        rows[-1]["base_protocol"] = [x.name for x in b.codes] if b else None
        print(f"{name} {r['eps0']:g}->{r['target']:g}: base {r['base']:.1f} {rows[-1]['base_protocol']} | proved with code {r['with_upper'] and round(r['with_upper'],1)} {r['protocol_proved']} -> {r['saving_proved'] and round(r['saving_proved'],2)}x", flush=True)
    out[name] = {"threshold_rigorous": c.threshold_eps, "rows": rows, "baseline": [x.name for x in base]}
    print(f"   rigorous threshold {c.threshold_eps:.3g}", flush=True)
json.dump(out, open(sys.argv[1], "w"), indent=1)
