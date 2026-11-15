"""Cost of the single-species paired-cap codes (census_paired_w8.py + isd_paired_w8.py) against [[1746,270,7]]_7, in both
species scenarios of species_cost.py.  The paired-cap codes enter as leading-order surrogates c eps^8 (L_8 from the census:
ESTIMATE ONLY, no proved tail yet); [[1746,270,7]]_7 with its proved weight-8/9 curve (species_w89.py) and its estimate.
usage: PYTHONPATH=.:scratchpad python scratchpad/paired_cap_cost.py L8_596 L8_606 OUT.json"""
import sys, json, math
L596, L606, OUT = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
sys.argv = ["species_cost.py"]
src = open("scratchpad/species_cost.py").read().replace("if __name__ == \"__main__\":\n    main()", "")
g = {"__name__": "paired_cap_cost"}
exec(compile(src, "species_cost.py", "exec"), g)
from dataclasses import dataclass
from qmsd.objective import LeadingOrderCurve, standard_library, DEFAULT_GRID
from qmsd.protocol import best_protocol
from w8_proved import upgrade


@dataclass
class Surr(LeadingOrderCurve):
    def f_surrogate(self, eps):
        return self.f(eps)

    @property
    def gamma(self):
        return math.log(self.n / self.k) / math.log(self.d)


new = {"[[596,90,8]]_7": Surr("[[596,90,8]]_7", 7, 596, 90, 8, L596, "paired cap |Y|=45, U only"),
       "[[606,80,8]]_7": Surr("[[606,80,8]]_7", 7, 606, 80, 8, L606, "paired cap |Y|=40, U only")}
sphere = upgrade(g["from_census"]("runs/weighted/p7_m4_sphere_1746_270_7_census.json"), 402)[0]
Tagged = g["Tagged"]
cands = {sphere.name: Tagged(sphere, {0: 672 - sphere.k, 1: 672, 2: 672})}
for nm, c in new.items(): cands[nm] = Tagged(c, {0: c.n})
tagged = g["tag_library"](); base_raw = standard_library(7); rows = []
for eps0, target in DEFAULT_GRID:
    ref = best_protocol(base_raw, eps0, target, max_depth=5); row = {"eps0": eps0, "target": target, "best_protocol_base": ref.cost}
    for sc in ("multi", "single"):
        b = g["answer"](g["search"](tagged, eps0, target, sc), target); row[sc] = {"base": b}
        line = f"{eps0:g}->{target:g} {sc:6s} base {b['cost']:8.2f} |"
        for nm, t in cands.items():
            est = g["answer"](g["search"](tagged + [t], eps0, target, sc, mode="surrogate"), target)
            ent = {"estimate": est, "saving_estimate": b["cost"] / est["cost"]}
            if nm == sphere.name:
                up = g["answer"](g["search"](tagged + [t], eps0, target, sc, mode="upper"), target)
                ent.update(proved=up, saving_proved=b["cost"] / up["cost"])
            row[sc][nm] = ent
            line += f" {nm} est {b['cost'] / est['cost']:.2f}x" + (f" (proved {ent['saving_proved']:.2f}x)" if "proved" in ent else "") + " |"
        print(line, flush=True)
    rows.append(row)
json.dump({"L8": {"[[596,90,8]]_7": L596, "[[606,80,8]]_7": L606}, "note": "paired-cap codes: leading-order estimate only",
           "rows": rows}, open(OUT, "w"), indent=1, default=str)
