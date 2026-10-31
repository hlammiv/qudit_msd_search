"""Rigorous protocol savings for the large weighted-support p = 7 codes (PERFORMANCE.md).
For each complete weight-d census, build the two-sided RigorousBoundCurve and compare the cheapest
protocol with the standard p = 7 library alone (exact curves) against the cheapest protocol that
may use the candidate with its PESSIMISTIC curve (proved saving) and with its OPTIMISTIC curve
(a bound on how much any proof could show).  Writes runs/weighted/rigorous_savings.json."""
import json, math
from qmsd.rigorous_curve import from_census, rigorous_saving, RigorousBoundCurve
from qmsd.objective import standard_library, LeadingOrderCurve, marginal_saving
base = standard_library(7)
out = {}
for path in ["runs/weighted/p7_m4_unit_2093_308_5_census.json", "runs/weighted/p7_m4_sphere_1746_270_6_census.json", "runs/weighted/p7_m4_sphere_1756_260_7_census.json"]:
    c = from_census(path); sur = LeadingOrderCurve(c.name, c.p, c.n, c.k, c.d, c.L_d)
    res = rigorous_saving(c, base)
    res.update({"A_d": c.A_d, "L_d": c.L_d, "c_eps": c.c_eps, "threshold_rigorous": c.threshold_eps,
                "f_hi_over_surrogate_at_1e-6": float(c.f_hi(1e-6) / c.f_surrogate(1e-6)),
                "f_hi_at_1e-6": float(c.f_hi(1e-6)), "surrogate_at_1e-6": float(c.f_surrogate(1e-6)),
                "surrogate_saving": marginal_saving(sur, base)})
    out[c.name] = res
    print(c.name, f"A_d={c.A_d} L_d={c.L_d:.3f} c_eps={c.c_eps:.3g} thr_rig={c.threshold_eps:.3g} f_hi/sur@1e-6={res['f_hi_over_surrogate_at_1e-6']:.3g}")
    for row in res["rows"]:
        print(f"   {row['eps0']:g} -> {row['target']:g}: base {row['base']:.1f}  proved {row['with_upper'] and round(row['with_upper'],1)} ({row['saving_proved'] and round(row['saving_proved'],2)}x)  optimistic {row['with_lower'] and round(row['with_lower'],1)} ({row['saving_optimistic'] and round(row['saving_optimistic'],2)}x)  uses={row['uses_proved']} {row['protocol_proved']} chain={row['eps_chain_proved'] and ['%.2e' % x for x in row['eps_chain_proved']]} recheck={row['reaches_target_check']}")
json.dump(out, open("runs/weighted/rigorous_savings.json", "w"), indent=1)
