# The corrected m=2 crux: the 2D-codeword cap (verified at the extremal arc)

2026-06-30 — "the real run" at the corrected crux (`D_PROOF_MAP.md` §4), **with a mid-run self-correction.**
Net outcome: the p=11 m=2 no-go is **verified at the extremal configuration** (true `d=4`, no crossing), the
binding object is confirmed to be the **2D full-span codeword** (the proof map's original crux — restored),
and a line/arc-counting shortcut I attempted **was refuted** by the classical arc value. The *general*
rigorous proof of the +1 gap remains **open**, with a now-sharpened, verified target.

## 1. Restatement (unchanged)
For p≤13 m=2, `G0 = shorten(RM_p(r_max,2), S)` (auto-triorthogonal; the moment hypothesis is inert),
`|S|=k`, `n_c=p²−k`. The scoped no-go ⟺ `d ≤ n_c/k` (γ≥1) ⟺ `k(d+1) ≤ p²`.

## 2. The attempted line+arc reduction — and why it is WRONG

I first reduced the crux to a pure incidence bound: `d ≤ d_RM − max_line` (**line bound L6, proved**) plus
`max_line ≥ 4` for every k=21 set (i.e. `m_3(AG(2,11)) ≤ 20`, the max ≤3-per-line set). My greedy/annealer/
ILP/SAT all capped construction at 20, which *looked* like `m_3=20`.

**This was a search artifact.** The literature is definitive: **`m_3(2,11) = 21`** (Marcugini,
*Maximal (n,3)-arcs in PG(2,11)*; Coolsaet–Sticker, *The complete (k,3)-arcs of PG(2,q), q≤13*), with only
**two** inequivalent (21,3)-arcs — rigid structures random search essentially never hits. Arc counting
forces any (21,3)-arc to miss `t_0 = 91 − t_3 ≥ 21` lines (`t_3≤70`), so it embeds affinely:
**`m_3(AG(2,11)) = 21`.** Hence `max_line = 3` **is** achievable at k=21, `d_lines = d_RM − 3 = 5`, and the
line bound alone does **not** close p=11. The line+arc route is refuted.

## 3. Verified verdict at the extremal arc (the decisive check)

HiGHS constructed an explicit (21,3)-arc (max_line=3, full rank, dim G0=7):
`S = [2,6,7,12,33,38,41,44,48,51,56,57,59,67,73,81,85,87,99,109,120]` (0-indexed AG(2,11) columns).

- Line bound: `d_lines = 5` (would give γ = 0.97 < 1 if tight).
- **True distance: `d = min_dependent_columns(G0) = 4`** ⇒ **`[[100,21,4]]₁₁`, γ = 1.126 ≥ 1.**

So a **2D full-span codeword caps `d` at 4 < `d_lines = 5`** — **no γ<1 crossing**, the no-go holds at the
extremal configuration, and the binding object is the 2D codeword (the proof map's original crux, **not** the
line bound). This matches the prior-session ILP finding and the workflow's "no p=11 crossing" stress-test.

## 4. The correct open crux (= the proof map's L-crux, restored)

> **Crux (open).** For p≤13 m=2, every full-rank puncture set `S` admits a **2D full-span codeword** in
> `ker(G0)` of weight `≤ n_c/k`, forcing `d ≤ n_c/k` (γ≥1) — even when `S` is a `max_line=3` extremal arc
> where the line bound gives `d_lines > n_c/k`. Verified at the p=11 extremal (21,3)-arc (`d=4`); the general
> statement (the DGM α=1 second-weight-class word capping `d`) is unproven.

The line/arc-counting shortcut does **not** reach this — it is genuinely a statement about the second
low-weight class of `RM_p(rtilde,2)` after puncturing, not about `max_line`.

## 5. What stands

- **Proved:** the line bound L6 (`d ≤ d_RM − max_line`); the elementary `max_line ≥ 3` bounds (second-moment
  `(p+k)/(p+1)`, pencil `⌈(k−1)/(p+1)⌉+1`) — now seen to be *consistent* with `m_3=21`, not a route to the +1.
- **Verified (not just evidenced):** at the p=11 extremal `max_line=3` arc, `true d = 4` ⇒ no crossing.
- **Open:** the general 2D-codeword cap (`d ≤ n_c/k`) — the true crux, harder than an arc bound.

## 6. Status and next steps

**p=11 m=2 no-go: holds, and verified at the extremal configuration; rigorous general proof OPEN** (crux =
2D-codeword cap). The clean "reduce to a classical arc bound" hope is closed off (`m_3(AG(2,11))=21`).
To finish, prove the 2D full-span codeword cap directly via the GRM (DGM) weight hierarchy — the α=1
second-weight class of the punctured dual — showing it punctures to `≤ n_c/k` for every full-rank `S` at
p≤13. This is a coding-theoretic (weight-hierarchy) problem, not a finite-geometry arc bound.

*Meta:* checking the literature (option 2) first was the right move — it exposed that the arc-bound target
was false before any lenore compute was spent chasing it.
