"""Does grading help?  p = 7 standard library (exact curves) against the same library plus graded
variants of its small codes (exact grade statistics, qmsd.graded), with and without the large p = 7
codes (leading-order surrogates from their censuses).  usage: graded_cost.py OUT.json"""
import sys, json, time
from qmsd.objective import standard_library, LeadingOrderCurve
from qmsd.protocol import best_protocol, _rm_k1, CYCLIC
from qmsd.cyclic import construction_A, construction_B
from qmsd.doubling import paper_codes
from qmsd.graded import GradedStats, GradedCurve
p = 7; t0 = time.time()
mats = {"[[6,1,2]]_7 RS": _rm_k1(p, 1)}
for label, cons, n, S in CYCLIC.get(p, []):
    if n <= 30: mats[label] = (construction_A if cons == "A" else construction_B)(p, n, S)
for label, G0, H1, c in paper_codes(p): mats[label] = (G0, H1)
tmax = {"[[7,1,3]]_7 cyc-B": 1, "[[12,1,3]]_7 cyc-A": 1, "[[20,1,4]]_7 dbl": 1, "[[28,1,5]]_7 cyc-B": 2}
graded = []
for name, t in tmax.items():
    if name not in mats: print("missing", name); continue
    st = GradedStats(name, *mats[name], p, t)
    configs = [(1, 0, True), (1, 1, False)] + ([(2, 0, True), (2, 1, True), (2, 2, False)] if t >= 2 else [])
    for (tt, g, rec) in configs:
        c = GradedCurve(st, tt, g, rec); graded.append(c)
        print(f"built {c.name} ({time.time()-t0:.0f}s): f(1e-2)={float(c.f(1e-2)):.3g} pacc_eff(1e-2)={c.pacc(1e-2):.3f}", flush=True)
base = standard_library(p)
big = [LeadingOrderCurve("[[1745,271,7]]_7", 7, 1745, 271, 7, 193.2), LeadingOrderCurve("[[1746,270,6]]_7", 7, 1746, 270, 6, 1.333)]
grid = [(1e-2, 1e-10), (1e-2, 1e-15), (1e-2, 1e-30), (3e-3, 1e-15), (1e-3, 1e-15), (1e-3, 1e-30)]
out = {"rows": []}
for eps0, tgt in grid:
    row = {"eps0": eps0, "target": tgt}
    for lab, pool in [("library", base), ("library+graded", base + graded), ("library+big", base + big), ("library+big+graded", base + big + graded)]:
        b = best_protocol(pool, eps0, tgt, max_depth=6)
        row[lab] = {"cost": b.cost, "protocol": [c.name for c in b.codes]} if b else None
    row["saving_graded_alone"] = row["library"]["cost"] / row["library+graded"]["cost"]
    row["saving_graded_with_big"] = row["library+big"]["cost"] / row["library+big+graded"]["cost"]
    out["rows"].append(row)
    print(f"{eps0:g}->{tgt:g}: library {row['library']['cost']:.1f} -> +graded {row['library+graded']['cost']:.1f} ({row['saving_graded_alone']:.2f}x) {row['library+graded']['protocol']} | with big codes {row['library+big']['cost']:.1f} -> {row['library+big+graded']['cost']:.1f} ({row['saving_graded_with_big']:.2f}x) {row['library+big+graded']['protocol']}", flush=True)
json.dump(out, open(sys.argv[1], "w"), indent=1)
