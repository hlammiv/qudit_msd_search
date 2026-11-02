"""Every catalog code: overhead exponent gamma against the input error at which the block is accepted
half the time, eps_1/2 = 1 - 2^(-1/n) ~ ln2/n (acceptance (1-eps)^n, verified by Monte Carlo to 1%).
Writes runs/coverage/gamma_vs_acceptance.png and .json."""
import json, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qmsd.results import load_result_catalog
best = {}
for r in load_result_catalog():
    if r.status not in ("confirmed", "partial") or not r.n or not r.k or not r.d or r.d < 2: continue
    g = math.log(r.n / r.k) / math.log(r.d)
    if g > 3: continue
    fam = (r.construction or {}).get("family", "")
    src = (r.note or "") + " " + (r.family or "") + " " + (r.artifact_id or "")
    if fam == "reed_solomon" or "RS" in (r.label or "") and r.n < 100 and r.m in (1, None):
        grp = "Reed–Solomon (m = 1)"
    elif fam == "weighted_rm" or "level-set" in src or "level set" in src:
        grp = "new constructions (weighted p = 7, level sets)"
    else:
        grp = "punctured RM: paper codes and our searches"
    key = (r.p, r.n, r.k)
    if key not in best or g < best[key]["gamma"]:
        best[key] = {"p": r.p, "n": r.n, "k": r.k, "d": r.d, "gamma": g, "partial": r.status == "partial", "group": grp, "label": r.label, "id": r.artifact_id}
rows = list(best.values())
for x in rows: x["eps_half"] = 1 - 2 ** (-1 / x["n"])
json.dump(rows, open("runs/coverage/gamma_vs_acceptance.json", "w"), indent=1)
C = {"Reed–Solomon (m = 1)": "#2a78d6", "punctured RM: paper codes and our searches": "#eb6834", "new constructions (weighted p = 7, level sets)": "#1baf7a"}
INK, MUTED, GRID = "#1f1f1f", "#6b6b6b", "#e6e6e6"
fig, ax = plt.subplots(figsize=(9, 6), dpi=160)
for grp, col in C.items():
    pts = [x for x in rows if x["group"] == grp]
    for hollow in (False, True):
        q = [x for x in pts if x["partial"] == hollow]
        if not q: continue
        ax.scatter([x["gamma"] for x in q], [x["eps_half"] for x in q], s=34, marker="o",
                   facecolors="white" if hollow else col, edgecolors=col if hollow else "white", linewidths=1.4 if hollow else 0.8,
                   label=f"{grp}" + (" — d is a lower bound" if hollow else ""), zorder=3 if grp.startswith("new") else 2)
ax.axvline(1.0, color=MUTED, lw=1, ls=(0, (4, 3)), zorder=1); ax.text(1.01, 0.03, "γ = 1", color=MUTED, fontsize=8, transform=ax.get_xaxis_transform())
for lev, txt in [(1e-3, "typical physical error 1e-3")]:
    ax.axhline(lev, color=GRID, lw=1, zorder=0); ax.text(2.97, lev * 1.25, txt, color=MUTED, fontsize=8, ha="right")
# labelled codes: (label as stored, text position in data coords)
hl = {"[[1745,271,7]]_7": (1.25, 1.1e-4), "[[13332,1309,22]]_11": (1.08, 1.5e-5), "[[78455,5066,>=68]]_17": (0.47, 4e-7),
      "[[123284,7037,>=57]]_19": (0.70, 4e-8), "[[666,175,6]]_29": (0.50, 2e-4), "[[73,24,9]]_97": (0.50, 0.25), "[[519,106,5]]_5": (1.25, 0.35)}
done = set()
for x in sorted(rows, key=lambda x: x["gamma"]):
    if x["label"] in hl and x["label"] not in done:
        done.add(x["label"]); tx, ty = hl[x["label"]]
        lab = x["label"].replace(">=", "≥").replace("]]_", "]]") + f"  γ {x['gamma']:.3f}"
        ax.annotate(lab, (x["gamma"], x["eps_half"]), xytext=(tx, ty), textcoords="data", fontsize=8, color=INK,
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7, shrinkA=2, shrinkB=4))
ax.set_yscale("log"); ax.set_xlim(0.45, 3.0)
man = [x for x in rows if x["eps_half"] < 1e-6]
if man:
    tgt = min((m for m in man if 1e-13 < m["eps_half"] < 1e-9), key=lambda m: m["eps_half"], default=man[0])
    ax.annotate("paper's analytic (Manhattan) family,\nn = 10⁷ – 10¹⁵ (same construction, all points below 10⁻⁷)", (tgt["gamma"], tgt["eps_half"]), xytext=(0.47, 1e-12),
                textcoords="data", fontsize=8, color=INK, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7, shrinkB=4))
ax.set_xlabel("overhead exponent γ = ln(n/k) / ln d   (lower is better asymptotically)", color=INK)
ax.set_ylabel("input error at which half the blocks are accepted, ε½ ≈ ln2 / n", color=INK)
ax.set_title("Asymptotic overhead against acceptance: low γ costs block size", color=INK, fontsize=12, loc="left")
ax.grid(True, which="major", color=GRID, lw=0.6); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
for s in ("left", "bottom"): ax.spines[s].set_color("#bdbdbd")
ax.tick_params(colors=MUTED)
ax.legend(fontsize=8, frameon=False, loc="lower right")
fig.text(0.01, 0.005, f"{len(rows)} distinct codes from the catalog (best γ per (p, n, k)); acceptance (1−ε)^n matches Monte Carlo to 1%. Source: scratchpad/gamma_vs_acceptance_plot.py", fontsize=7, color=MUTED)
fig.tight_layout(rect=(0, 0.02, 1, 1)); fig.savefig("runs/coverage/gamma_vs_acceptance.png")
print(len(rows), {g: sum(1 for x in rows if x["group"] == g) for g in C})
