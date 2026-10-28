"""Figures and reference numbers for the Monte Carlo guide: exact curves (qmsd.performance) with
Monte Carlo overlays (qmsd.mc_distill) for three small catalog codes."""
import json, math, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qmsd.performance import canonical_logical_rows, leading_order, profile_exact, record_matrices, threshold
from qmsd.reedmuller import r_max, rm_generator
from qmsd.results import load_result_catalog
from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.mc_distill import simulate
def rm_k1(p, m):
    r = r_max(m, p); G = rm_generator(r, m, p)
    return build_triorthogonal_code(p, m, r, [1], G=G)["G0"], canonical_logical_rows(G, [0], p)
recs = {r.artifact_id: r for r in load_result_catalog()}
r44 = next(r for r in load_result_catalog() if r.p == 7 and r.n == 44 and r.k == 5 and r.puncture_columns)
r18 = next(r for r in load_result_catalog() if r.p == 23 and r.n == 18 and r.k == 5)
codes = {"[[8,1,2]]_3": (3, rm_k1(3, 2)), "[[44,5,3]]_7": (7, record_matrices(r44)), "[[18,5,4]]_23": (23, record_matrices(r18))}
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "ink": "#2b2b2b", "muted": "#8a8a85", "grid": "#e6e6e3", "surface": "#fcfcfb"}
cols = {"[[8,1,2]]_3": C["blue"], "[[44,5,3]]_7": C["orange"], "[[18,5,4]]_23": C["aqua"]}
plt.rcParams.update({"font.size": 10, "figure.facecolor": C["surface"], "axes.facecolor": C["surface"], "axes.edgecolor": C["muted"]})
eps_grid = np.logspace(-3, math.log10(0.6), 120); mc_eps = [0.2, 0.1, 0.05, 0.02]; shots = 400_000
numbers = {}
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"wspace": 0.28})
for a in (a1, a2):
    a.grid(True, color=C["grid"]); a.set_axisbelow(True)
    for s in ("top", "right"): a.spines[s].set_visible(False)
for name, (p, (G0, H1)) in codes.items():
    prof = profile_exact(G0, H1, p, kinds=("avg",)); lo = leading_order(prof, "avg"); th = threshold(prof, "avg")
    ex = np.array([float(prof.eps_out(e, "avg")) for e in eps_grid]); pa = np.array([float(prof.p_accept(e)) for e in eps_grid])
    a1.plot(eps_grid, ex, color=cols[name], linewidth=2, label=f"{name}: exact")
    a1.plot(eps_grid, lo.c_eps * eps_grid ** lo.d, color=cols[name], linewidth=1, linestyle="--")
    a2.plot(eps_grid, pa, color=cols[name], linewidth=2, label=name)
    mcs = []
    for e in mc_eps:
        mc = simulate(G0, H1, p, e, shots=shots, seed=1); mcs.append((e, mc))
        a1.errorbar([e], [mc["eps_out_avg"]], yerr=[mc["eps_out_avg_err"]], fmt="o", color=cols[name], markersize=6, markeredgecolor=C["surface"], capsize=3)
        a2.errorbar([e], [mc["P_acc"]], yerr=[mc["P_acc_err"]], fmt="o", color=cols[name], markersize=6, markeredgecolor=C["surface"], capsize=3)
    if th.eps is not None:
        a1.plot([th.eps], [th.eps], "s", color=cols[name], markersize=8, markeredgecolor=C["ink"])
    numbers[name] = {"p": p, "n": int(G0.shape[1]), "k": int(H1.shape[0]), "dim_A": int(G0.shape[0]), "d": lo.d, "c_eps": float(lo.c_eps), "L_d": str(lo.L_d) if hasattr(lo, "L_d") else None,
                     "threshold_eps": float(th.eps) if th.eps is not None else None,
                     "exact": {str(e): {"P_acc": float(prof.p_accept(e)), "eps_out_avg": float(prof.eps_out(e, "avg"))} for e in (0.2, 0.1, 0.05, 0.02, 0.01)},
                     "mc": {str(e): {kk: (float(v) if isinstance(v, float) else v) for kk, v in mc.items()} for e, mc in mcs}}
    print(name, numbers[name]["d"], numbers[name]["c_eps"], numbers[name]["threshold_eps"], flush=True)
a1.plot(eps_grid, eps_grid, color=C["ink"], linewidth=1, linestyle=":", label="ε_out = ε (no gain)")
a1.set_xscale("log"); a1.set_yscale("log"); a1.set_xlabel("input error ε"); a1.set_ylabel("output error ε_out (per output)"); a1.set_ylim(1e-11, 1)
a1.legend(fontsize=8, frameon=False, loc="lower right"); a1.set_title("ε_out: exact, leading order (dashed), threshold (■), MC (●)", fontsize=9.5, loc="left")
a2.set_xscale("log"); a2.set_xlabel("input error ε"); a2.set_ylabel("acceptance probability P_acc"); a2.set_ylim(0, 1.02); a2.legend(fontsize=8, frameon=False, loc="lower left")
a2.set_title("P_acc: exact, MC (●, 4e5 shots)", fontsize=9.5, loc="left")
fig.savefig("docs/aanha_mc_guide/fig_curves.png", dpi=160, bbox_inches="tight")
# figure 2: per-output vs block for the k = 5 codes, and the MC resolution limit
fig, (b1, b2) = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw={"wspace": 0.28})
for a in (b1, b2):
    a.grid(True, color=C["grid"]); a.set_axisbelow(True)
    for s in ("top", "right"): a.spines[s].set_visible(False)
for name in ("[[44,5,3]]_7", "[[18,5,4]]_23"):
    p, (G0, H1) = codes[name]; prof = profile_exact(G0, H1, p, kinds=("avg",))
    b1.plot(eps_grid, [float(prof.eps_out(e, "avg")) for e in eps_grid], color=cols[name], linewidth=2, label=f"{name} per-output average (exact)")
    xs = [0.3, 0.2, 0.1, 0.05]; mcs = [simulate(G0, H1, p, e, shots=shots, seed=5) for e in xs]
    b1.errorbar(xs, [m["eps_out_block"] for m in mcs], yerr=[m["eps_out_block_err"] for m in mcs], fmt="s", color=cols[name], markersize=6, markeredgecolor=C["ink"], capsize=3, label=f"{name} block, any output wrong (MC)")
    b1.errorbar(xs, [m["eps_out_avg"] for m in mcs], yerr=[m["eps_out_avg_err"] for m in mcs], fmt="o", color=cols[name], markersize=6, markeredgecolor=C["surface"], capsize=3)
b1.set_xscale("log"); b1.set_yscale("log"); b1.set_xlabel("ε"); b1.set_ylabel("ε_out"); b1.set_ylim(1e-11, 1); b1.legend(fontsize=8, frameon=False, loc="lower right")
b1.set_title("k = 5: block (squares) vs per-output (line, dots)", fontsize=9.5, loc="left")
p, (G0, H1) = codes["[[8,1,2]]_3"]; prof = profile_exact(G0, H1, p, kinds=("avg",))
for sh, mk in ((20_000, "^"), (400_000, "o"), (4_000_000, "s")):
    xs, ys, es = [], [], []
    for e in (0.3, 0.1, 0.03, 0.01, 0.003):
        mc = simulate(G0, H1, p, e, shots=sh, seed=3); xs.append(e); ys.append(mc["eps_out_avg"]); es.append(mc["eps_out_avg_err"])
    b2.errorbar(xs, ys, yerr=es, fmt=mk, color=C["blue"], markersize=6, markeredgecolor=C["surface"], capsize=3, alpha=0.5 + 0.25 * (mk != "^"), label=f"{sh:,} shots")
b2.plot(eps_grid, [float(prof.eps_out(e, "avg")) for e in eps_grid], color=C["ink"], linewidth=1.5, label="exact")
b2.set_xscale("log"); b2.set_yscale("log"); b2.set_xlabel("ε"); b2.set_ylabel("ε_out, [[8,1,2]]_3"); b2.set_ylim(1e-7, 1); b2.legend(fontsize=8, frameon=False, loc="lower right")
b2.set_title("MC resolution: error ∝ 1/√(accepted · ε_out)", fontsize=9.5, loc="left")
fig.savefig("docs/aanha_mc_guide/fig_block_and_shots.png", dpi=160, bbox_inches="tight")
json.dump(numbers, open("docs/aanha_mc_guide/numbers.json", "w"), indent=1); print("done")
