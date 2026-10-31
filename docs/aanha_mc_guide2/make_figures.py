"""Figures and reference numbers for docs/aanha_mc_guide2 (thresholds of the large codes by
sampling, the two-round correlation test, biased noise).  Run from the repo root:
  PYTHONPATH=. OMP_NUM_THREADS=1 python docs/aanha_mc_guide2/make_figures.py"""
import json, math, itertools, time, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qmsd.results import load_result_catalog
from qmsd.performance import record_matrices, profile_exact
from qmsd.rigorous_curve import from_census
from qmsd.cyclic import construction_B
OUT = os.path.dirname(os.path.abspath(__file__)); NUM = {}
recs = {r.artifact_id: r for r in load_result_catalog()}

# ---------------- Task A: acceptance of the large codes and the window the sampler can reach
LARGE = [("[[1756,260,7]]_7", "imported-p7-n1756-k260-9d263f2240e4"), ("[[1746,270,7]]_7", "imported-p7-n1746-k270-c93a88f0735a"),
         ("[[666,175,6]]_29", "imported-p29-n666-k175-4f8ad426d19b"), ("[[766,195,>=6]]_31", "imported-p31-n766-k195-8cc6fbf9e4da"),
         ("[[424,105,6]]_23", "imported-p23-n424-k105-9d5204d5ccf4"), ("[[296,65,6]]_19", "imported-p19-n296-k65-211e1fbb4b1a")]
eps_grid = [1e-3, 2e-3, 3e-3, 5e-3, 7e-3, 1e-2]
NUM["acceptance"] = {}
fig, ax = plt.subplots(figsize=(6.2, 4))
es = np.logspace(-3.3, -1.7, 200)
for lab, aid in LARGE:
    n = recs[aid].n
    ax.plot(es, (1 - es) ** n, label=f"{lab} (n={n})")
    NUM["acceptance"][lab] = {f"{e:g}": (1 - e) ** n for e in eps_grid}
    NUM["acceptance"][lab]["shots_for_1000_accepted_at_3e-3"] = 1000 / (1 - 3e-3) ** n
for lvl, txt in [(1e-3, "1e-3: 1000 accepted per 1e6 shots"), (1e-5, "1e-5: 1000 accepted per 1e8 shots")]:
    ax.axhline(lvl, color="gray", ls=":", lw=0.8); ax.text(5.5e-4, lvl * 1.3, txt, fontsize=7, color="gray")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-8, 1.2); ax.set_xlabel(r"input error $\varepsilon$"); ax.set_ylabel(r"acceptance $(1-\varepsilon)^n$")
ax.set_title("Where postselected sampling is affordable"); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(f"{OUT}/fig_acceptance.png", dpi=160); plt.close(fig)

# ---------------- Task A: the rigorous curve, the surrogate and the diagonal for [[1756,260,7]]_7
c = from_census("runs/weighted/p7_m4_sphere_1756_260_7_census.json")
NUM["rigorous_1756"] = {"threshold_rigorous": c.threshold_eps, "c_eps": c.c_eps, "A_7": c.A_d, "L_7": c.L_d}
fig, ax = plt.subplots(figsize=(6.2, 4))
es = np.logspace(-4, -1.5, 120)
ax.plot(es, [float(c.f_hi(e)) for e in es], label=r"rigorous upper bound $f_{\rm hi}$ (Singleton tail)")
ax.plot(es, [float(c.f_surrogate(e)) for e in es], label=r"surrogate $c_\varepsilon\,\varepsilon^7$ (weight-7 term only)")
ax.plot(es, es, "k--", lw=0.8, label=r"$\varepsilon_{\rm out}=\varepsilon$")
ax.axvline(c.threshold_eps, color="C0", ls=":", lw=0.8); ax.text(c.threshold_eps * 1.05, 1e-20, "rigorous\nthreshold", fontsize=7, color="C0")
ax.axvspan(2e-3, 1.2e-2, color="orange", alpha=0.15); ax.text(2.2e-3, 1e-28, "sampler window\n(acceptance 1e-8 .. 3e-2)", fontsize=7, color="darkorange")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-30, 1); ax.set_xlabel(r"input error $\varepsilon$"); ax.set_ylabel(r"per-output error $\varepsilon_{\rm out}$")
ax.set_title(r"[[1756,260,7]]$_7$: what is known, what the sampler can add"); ax.legend(fontsize=7, loc="lower right"); fig.tight_layout(); fig.savefig(f"{OUT}/fig_rigorous_window.png", dpi=160); plt.close(fig)

# ---------------- Task B: how one accepted weight-7 error spreads over the outputs (census)
cen = json.load(open("runs/weighted/p7_m4_sphere_1756_260_7_census.json")); po = np.array(cen["per_output"])
NUM["census_spread_1756"] = {"outputs_hit_by_any_weight7_word": int((po > 0).sum()), "k": int(len(po)), "max_per_output": int(po.max()), "A_7": int(cen["A_7"]), "L_7": float(cen["L_d"])}
fig, axs = plt.subplots(1, 2, figsize=(8.5, 3.4))
axs[0].bar(range(len(po)), po, width=1.0, color="C3"); axs[0].set_xlabel("output index i (of 260)"); axs[0].set_ylabel("weight-7 words with $h_i\\cdot e\\neq0$")
axs[0].set_title("Each weight-7 word flips 6-12 outputs at once")
# routing schematic
ax = axs[1]; ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis("off")
for b in range(3):
    ax.add_patch(plt.Rectangle((0.3, 4.3 - 1.6 * b), 2.4, 1.2, fc="C0", alpha=0.25, ec="C0")); ax.text(1.5, 4.9 - 1.6 * b, f"round-1 block {b+1}\n(k outputs)", ha="center", va="center", fontsize=7)
ax.add_patch(plt.Rectangle((6.5, 3.9), 3.2, 1.6, fc="C2", alpha=0.25, ec="C2")); ax.text(8.1, 4.7, "round-2 block\ninputs from ONE block\n(correlated)", ha="center", va="center", fontsize=7)
ax.add_patch(plt.Rectangle((6.5, 0.6), 3.2, 1.6, fc="C1", alpha=0.25, ec="C1")); ax.text(8.1, 1.4, "round-2 block\ninputs shuffled across\nblocks (independent)", ha="center", va="center", fontsize=7)
for y0 in (4.6, 4.9, 5.2): ax.annotate("", xy=(6.5, 4.7), xytext=(2.7, y0), arrowprops=dict(arrowstyle="->", color="C2", lw=0.8))
for y0 in (4.9, 3.3, 1.7): ax.annotate("", xy=(6.5, 1.4), xytext=(2.7, y0), arrowprops=dict(arrowstyle="->", color="C1", lw=0.8))
ax.set_title("The two routings of exercise B3"); fig.tight_layout(); fig.savefig(f"{OUT}/fig_correlation.png", dpi=160); plt.close(fig)

# ---------------- Task B reference run: pair-correlation ratio on [[111,10,6]]_11 (own sampler, per-shot indicators)
def sample_failures(G0, H1, p, eps, shots, seed=0, batch=100_000, weights=None):
    """Return the (accepted, k) boolean failure matrix and the number of shots drawn."""
    G0 = np.asarray(G0) % p; H1 = np.asarray(H1) % p; n = G0.shape[1]; rng = np.random.default_rng(seed); out = []
    vals = np.arange(1, p)
    for s in range(0, shots, batch):
        b = min(batch, shots - s)
        v = rng.choice(vals, size=(b, n), p=weights) if weights is not None else rng.integers(1, p, size=(b, n))
        j = np.where(rng.random((b, n)) < eps, v, 0)
        ok = ~np.any((j @ G0.T) % p, axis=1)
        out.append(((j[ok] @ H1.T) % p) != 0)
    return np.concatenate(out), shots
# a d = 6 code ([[111,10,6]]_11) gives no failure at all in 3e5 accepted shots at eps <= 0.03; the test needs
# a code whose failures are frequent: the d = 3, k = 22 code [[99,22,3]]_11 and the paper's [[200,43,3]]_3
NUM["pair_correlation"] = {}
for lab, aid, p_, epss, shots in [("[[99,22,3]]_11", "imported-p11-n99-k22-1a07d37d9bc4", 11, (0.01, 0.03), 2_000_000),
                                  ("[[200,43,3]]_3", "paper-table3-p3-m5-n200-k43", 3, (0.005, 0.01), 2_000_000)]:
  G0, H1 = record_matrices(recs[aid]); NUM["pair_correlation"][lab] = {}
  for eps in epss:
    t0 = time.time(); F, shots = sample_failures(G0, H1, p_, eps, shots, seed=3)
    m = F.mean(axis=0); pair = (F[:, :, None] & F[:, None, :]).mean(axis=0)
    off = ~np.eye(F.shape[1], dtype=bool); ratio = pair[off] / np.outer(m, m)[off]
    nf = F.sum(axis=1); kk = F.shape[1]
    NUM["pair_correlation"][lab][f"{eps:g}"] = {"accepted": int(F.shape[0]), "P_acc": F.shape[0] / shots, "eps_out_avg": float(m.mean()), "eps_out_block": float(F.any(axis=1).mean()),
        "pair_ratio_mean": float(np.nanmean(ratio)), "pair_ratio_min": float(np.nanmin(ratio)), "pair_ratio_max": float(np.nanmax(ratio)),
        "mean_failures_per_failing_block": float(nf[nf > 0].mean()) if (nf > 0).any() else None,
        "P_two_or_more_fail_given_one_fails": float((nf >= 2).sum() / max((nf >= 1).sum(), 1)),
        "independent_prediction_P_two_or_more_given_one": float(1 - (1 - m.mean()) ** (kk - 1)), "seconds": round(time.time() - t0)}
    print("pair corr", lab, eps, NUM["pair_correlation"][lab][f"{eps:g}"], flush=True)

# ---------------- Task C: biased noise, exact by enumeration, [[7,1,3]]_7 cyc-B and [[6,1,?]]_7 RS
def exact_biased(G0, H1, p, eps, weights):
    """Exact per-output error and acceptance under P(Z^j) = eps * weights[j-1] by full enumeration (p^n vectors)."""
    G0 = np.asarray(G0) % p; H1 = np.asarray(H1) % p; n = G0.shape[1]
    allv = np.array(list(itertools.product(range(p), repeat=n)), dtype=np.int64)
    ok = ~np.any((allv @ G0.T) % p, axis=1); acc = allv[ok]
    w = np.array([1 - eps] + [eps * x for x in weights]); pr = np.prod(w[acc], axis=1)
    fail = ((acc @ H1.T) % p != 0).any(axis=1)
    return float(pr.sum()), float(pr[fail].sum() / pr.sum())
def bias_weights(p, b):
    """One nonzero power carries b + (1-b)/(p-1), the rest (1-b)/(p-1) each; b = 0 is the twirled model."""
    base = (1 - b) / (p - 1); return [b + base] + [base] * (p - 2)
NUM["biased"] = {}
G7, H7 = construction_B(7, 8, {1, 4, 7}); print("[[7,1,3]]_7 shapes", np.asarray(G7).shape, np.asarray(H7).shape)
from qmsd.rs_family import r_max_rs
from qmsd.reedmuller import rm_generator
from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.performance import canonical_logical_rows
G = rm_generator(r_max_rs(7), 1, 7); G6 = build_triorthogonal_code(7, 1, r_max_rs(7), [1], G=G)["G0"]; H6 = canonical_logical_rows(G, [0], 7)
bs = [0, 0.25, 0.5, 0.75, 0.9, 1.0]
fig, ax = plt.subplots(figsize=(6.2, 3.8))
for lab, (g, h) in {"[[7,1,3]]_7 cyc-B": (G7, H7), "[[6,1,?]]_7 RS": (G6, H6)}.items():
    NUM["biased"][lab] = {}
    for eps in (1e-2, 3e-2):
        ys = []
        for b in bs:
            pa, eo = exact_biased(g, h, 7, eps, bias_weights(7, b)); ys.append(eo); NUM["biased"][lab][f"eps={eps:g},b={b:g}"] = {"P_acc": pa, "eps_out": eo}
        ax.plot(bs, ys, "o-", label=f"{lab}, $\\varepsilon$={eps:g}")
ax.set_yscale("log"); ax.set_xlabel("bias b (0 = twirled, 1 = a single Z power)"); ax.set_ylabel(r"$\varepsilon_{\rm out}$ (exact, enumeration)")
ax.set_title("Biased noise changes the output error of small codes"); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(f"{OUT}/fig_biased.png", dpi=160); plt.close(fig)
# which power is favoured matters: scan the favoured power for [[7,1,3]]_7 at b = 0.9
NUM["biased_which_power"] = {}
for jstar in range(1, 7):
    wts = [(1 - 0.9) / 6] * 6; wts[jstar - 1] += 0.9
    NUM["biased_which_power"][f"j*={jstar}"] = exact_biased(G7, H7, 7, 1e-2, wts)[1]
json.dump(NUM, open(f"{OUT}/numbers.json", "w"), indent=1); print(json.dumps({k: v for k, v in NUM.items() if k != "acceptance"}, indent=1)[:3000])
