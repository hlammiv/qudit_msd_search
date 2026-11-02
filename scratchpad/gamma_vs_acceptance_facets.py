"""Small multiples of gamma vs eps_1/2 (input error at 50% acceptance), one panel per prime, so every
comparison is within one qudit dimension.  Data: runs/coverage/gamma_vs_acceptance.json (from
gamma_vs_acceptance_plot.py).  Writes runs/coverage/gamma_vs_acceptance_by_prime.png."""
import json, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
rows = json.load(open("runs/coverage/gamma_vs_acceptance.json"))
C = {"Reed–Solomon (m = 1)": "#2a78d6", "punctured RM: paper codes and our searches": "#eb6834", "new constructions (weighted p = 7, level sets)": "#1baf7a"}
SHORT = {"Reed–Solomon (m = 1)": "Reed–Solomon", "punctured RM: paper codes and our searches": "punctured RM (paper + our searches)", "new constructions (weighted p = 7, level sets)": "new constructions (weighted, level sets)"}
INK, MUTED, GRID = "#1f1f1f", "#6b6b6b", "#e6e6e6"
primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
fig, axs = plt.subplots(2, 5, figsize=(16, 7.2), dpi=150, sharey=True)
for ax, p in zip(axs.flat, primes):
    pts = [x for x in rows if x["p"] == p and x["gamma"] <= 3]
    for grp, col in C.items():
        for hollow in (False, True):
            q = [x for x in pts if x["group"] == grp and x["partial"] == hollow]
            if q:
                ax.scatter([x["gamma"] for x in q], [x["eps_half"] for x in q], s=26, facecolors="white" if hollow else col,
                           edgecolors=col if hollow else "white", linewidths=1.2 if hollow else 0.6, zorder=3 if grp.startswith("new") else 2)
    small = [x for x in pts if x["n"] <= 10 ** 6]
    if small:
        b = min(small, key=lambda x: x["gamma"])
        ax.annotate(b["label"].replace(">=", "≥").replace("]]_", "]]") + f"\nγ {b['gamma']:.3f}  (best with n ≤ 10⁶)", (b["gamma"], b["eps_half"]),
                    xytext=(1.55, 3e-12), textcoords="data", fontsize=7, color=INK, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkB=3))
    ax.axvline(1.0, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.axhline(1e-3, color=GRID, lw=0.8, zorder=0)
    ax.set_yscale("log"); ax.set_xlim(0.45, 3.0); ax.set_title(f"p = {p}   ({len(pts)} codes)", fontsize=10, color=INK, loc="left")
    ax.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color("#bdbdbd")
    ax.tick_params(colors=MUTED, labelsize=8)
for ax in axs[1]: ax.set_xlabel("γ = ln(n/k)/ln d", color=INK, fontsize=9)
for ax in axs[:, 0]: ax.set_ylabel("ε at 50% acceptance (≈ ln2/n)", color=INK, fontsize=9)
handles = [plt.Line2D([], [], marker="o", ls="", markerfacecolor=c, markeredgecolor="white", markersize=7, label=SHORT[g]) for g, c in C.items()]
handles.append(plt.Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor=MUTED, markersize=7, label="distance is a lower bound"))
fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.94))
fig.suptitle("γ against acceptance within each prime: lower γ is bought with block size", y=0.99, fontsize=12, color=INK, x=0.01, ha="left")
fig.text(0.01, 0.005, "Best γ per (p, n, k) among confirmed/partial catalog codes; dashed line γ = 1; grey line ε = 1e-3. Source: scratchpad/gamma_vs_acceptance_facets.py", fontsize=7, color=MUTED)
fig.tight_layout(rect=(0, 0.02, 1, 0.90)); fig.savefig("runs/coverage/gamma_vs_acceptance_by_prime.png")
print("saved")
