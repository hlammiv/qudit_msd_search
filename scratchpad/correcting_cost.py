"""Protocol costs when the large codes may also be run in detect-and-correct mode (t = 0 .. floor((d-1)/2),
CorrectingBoundCurve, census-free rigorous bounds) against the standard library + RS family at the
same prime.  Every t is offered to the optimizer as a separate curve; it picks the best.
usage: correcting_cost.py OUT.json"""
import sys, json
from qmsd.objective import standard_library
from qmsd.protocol import CodeCurve, best_protocol
from qmsd.rs_family import rs_profile, rs_parameters, r_max_rs
from qmsd.rigorous_curve import CorrectingBoundCurve
codes = [("[[13332,1309,22]]_11", 11, 13332, 1309, 22, [0, 2, 4, 6, 8, 10]), ("[[78455,5066,68]]_17", 17, 78455, 5066, 68, [0, 5, 10, 20, 30]),
         ("[[123284,7037,57]]_19", 19, 123284, 7037, 57, [0, 5, 10, 20, 28]), ("[[1745,271,7]]_7", 7, 1745, 271, 7, [0, 1, 2, 3])]
grid = [(1e-2, 1e-15), (1e-2, 1e-30), (1e-3, 1e-15), (1e-3, 1e-30), (1e-3, 1e-60)]
out = {}
for name, p, n, k, d, ts in codes:
    base = standard_library(p)
    for kk in range(1, r_max_rs(p) + 1):
        nn, _, dd = rs_parameters(p, kk)
        if dd >= 2: base.append(CodeCurve(f"[[{nn},{kk},{dd}]]_{p} RS", p, nn, kk, rs_profile(p, kk)))
    curves = [CorrectingBoundCurve(f"{name} t={t}", p, n, k, d, t) for t in ts]
    only_det = [curves[0]]
    out[name] = []
    for eps0, tgt in grid:
        b = best_protocol(base, eps0, tgt, max_depth=6)
        w0 = best_protocol(base + only_det, eps0, tgt, max_depth=6, mode="upper")
        wt = best_protocol(base + curves, eps0, tgt, max_depth=6, mode="upper")
        row = {"eps0": eps0, "target": tgt, "base": b.cost, "detect_only": w0.cost, "with_correction": wt.cost,
               "saving_detect": b.cost / w0.cost, "saving_correct": b.cost / wt.cost, "protocol": [c.name for c in wt.codes], "eps_chain": [float(x) for x in wt.eps]}
        out[name].append(row)
        print(f"{name} {eps0:g}->{tgt:g}: base {b.cost:.1f} | detect-only {w0.cost:.1f} ({row['saving_detect']:.2f}x) | with correction {wt.cost:.1f} ({row['saving_correct']:.2f}x) {[c.name.split(' ')[-1] if 't=' in c.name else c.name.split(' ')[0] for c in wt.codes]}", flush=True)
    for c in curves: print(f"   {c.name}: rigorous threshold {c.threshold_eps}", flush=True)
json.dump(out, open(sys.argv[1], "w"), indent=1)
