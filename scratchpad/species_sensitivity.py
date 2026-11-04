"""How far can raw U^2, U^3 differ from raw U before the p = 7 sphere-charge savings go?  Raw U at eps0, cost 1;
raw U^2, U^3 at r*eps0, cost w.  Uses the search in species_cost.py.  Writes runs/weighted/species_sensitivity.json."""
import json, sys
sys.path.insert(0, "scratchpad")
import species_cost as sc
from qmsd.rigorous_curve import from_census

orig = sc.Pareto


def run(pool, eps0, target, r, w, mode):
    """Multi-species search with modified raw points: seed via a Pareto that rewrites the raw entries."""
    calls = {"n": 0}

    class Seeded(orig):
        def __init__(self):
            super().__init__()
            self.idx = calls["n"]
            calls["n"] += 1

        def add(self, eps, cost, how):
            if how == "raw" and self.idx in (1, 2):
                return super().add(eps * r, w, how)
            return super().add(eps, cost, how)

    sc.Pareto = Seeded
    try:
        return sc.answer(sc.search(pool, eps0, target, "multi", mode=mode), target)
    finally:
        sc.Pareto = orig


def main():
    tagged = sc.tag_library()
    c = from_census("runs/weighted/p7_m4_sphere_1746_270_7_census.json")
    t = sc.Tagged(c, {0: 402, 1: 672, 2: 672})
    u = from_census("runs/weighted/p7_m4_unit_2093_308_5_census.json")
    tu = sc.Tagged(u, {0: u.n})
    out = []
    for eps0, target in ((1e-2, 1e-30), (1e-3, 1e-30)):
        for r in (0.7, 1.0, 1.3, 2.0, 3.0):
            for w in (1, 2, 3, 5):
                b = run(tagged, eps0, target, r, w, "upper")
                row = {"eps0": eps0, "target": target, "r": r, "w": w, "base": b["cost"]}
                for name, ct in (("1746", t), ("2093", tu)):
                    up = run(tagged + [ct], eps0, target, r, w, "upper")
                    es = run(tagged + [ct], eps0, target, r, w, "surrogate")
                    row[name] = {"proved": b["cost"] / up["cost"], "estimate": b["cost"] / es["cost"]}
                out.append(row)
                print(f"{eps0:g} r={r} w={w} base {b['cost']:.1f}: 1746 proved {row['1746']['proved']:.2f} est {row['1746']['estimate']:.2f}; 2093 proved {row['2093']['proved']:.2f} est {row['2093']['estimate']:.2f}", flush=True)
    json.dump(out, open("runs/weighted/species_sensitivity.json", "w"), indent=1)


if __name__ == "__main__":
    main()
