"""gamma vs p and block size vs p: arXiv:2510.10852 (asymptotic optimum gamma_0(p), Table 2 Manhattan
codes, Table 3 search codes) against our catalog's best certified code per prime (COVERAGE.md)."""
import json, math, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
gamma0 = {3: .63215, 5: .55915, 7: .50786, 11: .44108, 13: .41785, 17: .38273, 19: .36901, 23: .34663, 29: .32151, 31: .31468, 37: .29744, 41: .28798,
          43: .28373, 47: .27600, 53: .26602, 59: .25752, 61: .25495, 67: .24793, 71: .24372, 73: .24174, 79: .23624, 83: .23290, 89: .22829, 97: .22278}
# paper Table 2 (Manhattan family, smallest gamma<1 of the construction): p -> (n, gamma)
paper_manh = {3: (1852445880782154, .990), 5: (152186472515, .980), 7: (96448935471, .990), 11: (18874416, .964), 13: (60848853, .975),
              17: (77540, .946), 19: (121470, .967), 23: (17, .948)}
# paper Table 3 (search codes, best gamma per p)
paper_search = {3: (206, 1.239), 5: (519, 0.987)}
# ours (COVERAGE.md): p -> (n, gamma, label, kind)
ours = {3: (206, 1.239, "[[206,37,4]]", "paper puncture set"), 5: (519, 0.987, "[[519,106,5]]", "paper code"),
        7: (1756, 0.982, "[[1756,260,7]]", "new"), 11: (100, 1.126, "[[100,21,4]]", "new (m=2 search)"), 13: (141, 1.166, "[[141,28,4]]", "new (m=2 search)"),
        17: (237, 0.847, "[[237,52,6]]", "new"), 19: (296, 0.846, "[[296,65,6]]", "new"),
        23: (429, .813, "[[429,100,6]]", "new (m=2 search)"), 29: (22, .826, "", "RS"), 31: (25, .887, "", "RS"), 37: (29, .800, "", "RS"), 41: (31, .703, "", "RS"), 43: (34, .742, "", "RS"),
        47: (36, .662, "", "RS"), 53: (40, .627, "", "RS"), 59: (45, .600, "", "RS"), 61: (47, .622, "", "RS"), 67: (51, .596, "", "RS"), 71: (53, .555, "", "RS"),
        73: (56, .573, "", "RS"), 79: (60, .553, "", "RS"), 83: (62, .521, "", "RS"), 89: (67, .507, "", "RS"), 97: (73, .506, "", "RS")}
C = {"yellow": "#eda100", "blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "violet": "#6250d6", "ink": "#2b2b2b", "muted": "#8a8a85", "grid": "#e6e6e3", "surface": "#fcfcfb"}
plt.rcParams.update({"font.size": 10, "axes.edgecolor": C["muted"], "axes.labelcolor": C["ink"], "xtick.color": C["ink"], "ytick.color": C["ink"], "figure.facecolor": C["surface"], "axes.facecolor": C["surface"]})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.5, 8.8), sharex=True, gridspec_kw={"height_ratios": [1.15, 1], "hspace": 0.12})
for a in (a1, a2):
    a.grid(True, color=C["grid"], linewidth=0.8); a.set_axisbelow(True)
    for s in ("top", "right"): a.spines[s].set_visible(False)
# --- panel 1: gamma vs p
a1.plot(primes, [gamma0[p] for p in primes], color=C["blue"], linewidth=2, label="paper: asymptotic optimum γ₀(p) (Fig. 6, n → ∞)")
a1.plot(list(paper_manh), [v[1] for v in paper_manh.values()], "o", color=C["orange"], markersize=8, markeredgecolor=C["surface"], markeredgewidth=1.5, label="paper: Manhattan codes (Table 2, n up to 10¹⁵)")
a1.plot(list(paper_search), [v[1] for v in paper_search.values()], "s", color=C["orange"], markersize=8, markeredgecolor=C["surface"], markeredgewidth=1.5, label="paper: search codes (Table 3)")
new = [p for p in primes if ours[p][3].startswith("new")]; rs = [p for p in primes if ours[p][3] == "RS"]; same = [p for p in primes if ours[p][3].startswith("paper")]
a1.plot(rs, [ours[p][1] for p in rs], "D", color=C["aqua"], markersize=7, markeredgecolor=C["surface"], markeredgewidth=1.5, label="our catalog: Reed–Solomon family, min-γ member (p ≥ 23)")
a1.plot(new, [ours[p][1] for p in new], "D", color=C["violet"], markersize=9, markeredgecolor=C["surface"], markeredgewidth=1.5, label="our catalog: new search codes")
a1.plot(same, [ours[p][1] for p in same], "D", color=C["muted"], markersize=7, markeredgecolor=C["surface"], markeredgewidth=1.5, label="our catalog = paper's code (p = 3, 5)")
k1best = {3: (37, 2.242, "[[37,1,5]]₃ doubling"), 5: (16, 2.000, "[[16,1,4]]₅ cyclic"), 7: (7, 1.771, "[[7,1,3]]₇ cyclic")}
a1.plot(list(k1best), [v[1] for v in k1best.values()], "^", color=C["yellow"], markersize=9, markeredgecolor=C["surface"], markeredgewidth=1.5, label="2026 papers: best k = 1 cyclic / doubling code (arXiv:2610.08012, 2610.29457)")
for p, v in k1best.items(): a1.annotate(v[2], (p, v[1]), textcoords="offset points", xytext=(9, -4), fontsize=7.5, color=C["ink"])
a2.plot(list(k1best), [v[0] for v in k1best.values()], "^", color=C["yellow"], markersize=9, markeredgecolor=C["surface"], markeredgewidth=1.5, label="2026 papers: best k = 1 code")
a1.axhline(1.0, color=C["ink"], linewidth=1, linestyle="--"); a1.text(97, 1.03, "γ = 1", ha="right", va="bottom", color=C["ink"], fontsize=9)
a1.text(0.985, 0.62, "our new codes\n" + "\n".join(f"p = {p}: {ours[p][2]}  γ = {ours[p][1]:.3f}" for p in new), transform=a1.transAxes, ha="right", va="top", fontsize=8, color=C["ink"], bbox=dict(boxstyle="round,pad=0.4", facecolor=C["surface"], edgecolor=C["grid"]))
a1.set_ylabel("γ = ln(n/k) / ln d"); a1.set_ylim(0.15, 2.45); a1.legend(loc="upper right", fontsize=8, frameon=False)
a1.set_title("Overhead exponent γ by qudit dimension p: arXiv:2510.10852 vs our catalog (best certified code per prime)", fontsize=10.5, loc="left", color=C["ink"])
# --- panel 2: block size vs p (log)
a2.plot(list(paper_manh), [v[0] for v in paper_manh.values()], "o", color=C["orange"], markersize=8, markeredgecolor=C["surface"], markeredgewidth=1.5, label="paper: Manhattan codes n_min (Fig. 2)")
a2.plot(list(paper_search), [v[0] for v in paper_search.values()], "s", color=C["orange"], markersize=8, markeredgecolor=C["surface"], markeredgewidth=1.5, label="paper: search codes")
a2.plot(rs, [ours[p][0] for p in rs], "D", color=C["aqua"], markersize=7, markeredgecolor=C["surface"], markeredgewidth=1.5, label="our catalog: RS family")
a2.plot(new, [ours[p][0] for p in new], "D", color=C["violet"], markersize=9, markeredgecolor=C["surface"], markeredgewidth=1.5, label="our catalog: new search codes")
a2.plot(same, [ours[p][0] for p in same], "D", color=C["muted"], markersize=7, markeredgecolor=C["surface"], markeredgewidth=1.5)
a2.set_yscale("log"); a2.set_ylabel("block size n (qudits)"); a2.set_xlabel("p"); a2.set_xlim(0, 100); a2.legend(loc="upper right", fontsize=8, frameon=False)
a2.set_title("Block size of the best γ < 1 code (γ > 1 shown for p = 3, 11, 13)", fontsize=10.5, loc="left", color=C["ink"])
fig.savefig("runs/coverage/gamma_vs_p.png", dpi=160, bbox_inches="tight"); print("saved")
