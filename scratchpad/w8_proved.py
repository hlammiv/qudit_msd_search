"""Proved savings with the rigorous weight-8 bound for the d = 7 sphere-charge codes.
N_8 <= 6 C(b, 8) [block words: an 8-set carries a subcode of dim <= 2 (Singleton, d >= 7), and a
[8,2,7] code is MDS with no weight-8 word, so <= p-1 full-support words per 8-set]
     + A_8(RM_7(20,4)) [mixed words with all 8 points off the block: u = (g-1)(g-6) e_N]
     + A_5 + A_6 + A_7 (RM_7(20,4)) [mixed words with 1-3 block points: the block completion is unique, d >= 7].
L_8 <= N_8.  Writes runs/weighted/w8_proved.json."""
import json
from qmsd.rigorous_curve import RigorousBoundCurve, from_census, rigorous_saving
from qmsd.objective import standard_library

from qmsd.rigorous_curve import sphere_weight8_bound, RM20_A8
assert RM20_A8 == json.load(open("runs/weighted/rm20_w8_count.json"))["A8_RM20"]

def upgrade(c, block):
    """the census curve c with the proved weight-8 bound (library: RigorousBoundCurve.heavier)"""
    assert c.d == 7
    N8 = sphere_weight8_bound(block)
    return RigorousBoundCurve(c.name, c.p, c.n, c.k, c.d, c.L_d, c.A_d, heavier={8: (N8, N8)}), N8

if __name__ == "__main__":
    base = standard_library(7); out = {}
    for path, block in (("runs/weighted/p7_m4_sphere_1746_270_7_census.json", 402),
                        ("runs/weighted/p7_m4_sphere_1756_260_7_census.json", 412)):
        c = from_census(path); w, N8 = upgrade(c, block)
        old = rigorous_saving(c, base)["rows"]; new = rigorous_saving(w, base)["rows"]
        out[c.name] = {"N8_bound": N8, "threshold_old": c.threshold_eps, "threshold_new": w.threshold_eps,
                       "f_hi_1e-6_old": float(c.f_hi(1e-6)), "f_hi_1e-6_new": float(w.f_hi(1e-6)),
                       "surrogate_1e-6": float(c.f_surrogate(1e-6)), "rows": []}
        print(c.name, f"N8 <= {N8:.3e}; rigorous threshold {c.threshold_eps:.4g} -> {w.threshold_eps:.4g}; f_hi(1e-6) {float(c.f_hi(1e-6)):.2e} -> {float(w.f_hi(1e-6)):.2e} (surrogate {float(c.f_surrogate(1e-6)):.2e})")
        for a, b in zip(old, new):
            out[c.name]["rows"].append({"eps0": a["eps0"], "target": a["target"], "proved_old": a["saving_proved"], "proved_new": b["saving_proved"],
                                        "optimistic": a["saving_optimistic"], "protocol_new": b["protocol_proved"], "chain_check": b["reaches_target_check"]})
            print(f"   {a['eps0']:g} -> {a['target']:g}: proved {a['saving_proved']:.3f} -> {b['saving_proved']:.3f} (optimistic {a['saving_optimistic']:.3f}) {b['protocol_proved']} check={b['reaches_target_check']}")
    json.dump(out, open("runs/weighted/w8_proved.json", "w"), indent=1)
