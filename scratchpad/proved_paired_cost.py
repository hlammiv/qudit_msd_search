"""Proved costs of the paired-cap codes from an EXACT weight-8 census (census_paired_w8.py with the w8_mitm44v2.py part (ii)):
RigorousBoundCurve(A_8, L_8) with the weight-9/10 bounds of paired_heavier.py and the generic tail of rigorous_curve.py, in both species scenarios of species_cost.py, next to the
proved sphere [[1746,270,7]]_7 curve (weight-8/9 bounds).  Also the proved threshold (largest eps with f_hi(eps) < eps).
usage: PYTHONPATH=.:scratchpad python scratchpad/proved_paired_cost.py OUT.json CENSUS1.json [CENSUS2.json ...]"""
import sys, json
OUT, CENS = sys.argv[1], sys.argv[2:]
sys.argv = ["species_cost.py"]
src = open("scratchpad/species_cost.py").read().replace("if __name__ == \"__main__\":\n    main()", "")
g = {"__name__": "proved_paired_cost"}
exec(compile(src, "species_cost.py", "exec"), g)
from qmsd.rigorous_curve import RigorousBoundCurve
from qmsd.objective import standard_library, DEFAULT_GRID
from w8_proved import upgrade
from paired_heavier import heavier_bounds

Tagged = g["Tagged"]
sphere = upgrade(g["from_census"]("runs/weighted/p7_m4_sphere_1746_270_7_census.json"), 402)[0]
cands = {sphere.name: Tagged(sphere, {0: 672 - sphere.k, 1: 672, 2: 672})}
info = {"thresholds": {sphere.name: sphere.threshold_eps}}
for path in CENS:
    c = json.load(open(path)); assert c["part_ii_complete"], f"{path}: census not exact"
    nm = f"[[{c['n']},{c['k']},{c['d']}]]_7"
    hv = heavier_bounds(343 - c["k"] // 2, c["dprime_words"])
    cur = RigorousBoundCurve(nm, 7, c["n"], c["k"], c["d"], float(c["L8_found"]), int(c["A8_found"]), heavier=hv,
                             note="paired cap, exact weight-8 census + proved weight-9/10 bounds + generic tail")
    info.setdefault("heavier", {})[nm] = {str(w): v[0] for w, v in hv.items()}
    cands[nm] = Tagged(cur, {0: cur.n})
    info["thresholds"][nm] = cur.threshold_eps
    print(f"{nm}: A_8 {cur.A_d}, L_8 {cur.L_d:.1f}, gamma {cur.gamma:.3f}, proved threshold {cur.threshold_eps:.4%}", flush=True)
print(f"{sphere.name}: proved threshold {sphere.threshold_eps:.4%}", flush=True)
tagged = g["tag_library"](); rows = []
for eps0, target in DEFAULT_GRID:
    row = {"eps0": eps0, "target": target}
    for sc in ("multi", "single"):
        b = g["answer"](g["search"](tagged, eps0, target, sc), target); row[sc] = {"base": b}
        line = f"{eps0:g}->{target:g} {sc:6s} base {b['cost']:8.2f} |"
        for nm, t in cands.items():
            up = g["answer"](g["search"](tagged + [t], eps0, target, sc, mode="upper"), target)
            est = g["answer"](g["search"](tagged + [t], eps0, target, sc, mode="surrogate"), target)
            row[sc][nm] = {"proved": up, "estimate": est, "saving_proved": b["cost"] / up["cost"], "saving_estimate": b["cost"] / est["cost"]}
            line += f" {nm} proved {b['cost'] / up['cost']:.2f}x est {b['cost'] / est['cost']:.2f}x |"
        print(line, flush=True)
    rows.append(row)
info["rows"] = rows
json.dump(info, open(OUT, "w"), indent=1, default=str)
