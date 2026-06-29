# Which qudit dimension *p* minimizes magic-state-distillation overhead?

A synthesis over the codes and cost models produced in this project (Tier-5 analysis).
All numbers below are reproduced from `rs_family.py`, `qmsd/asymptotics.py`,
`qmsd/distillation.py`, and the closed forms verified against `qmsd.weightdist` in this
writeup. No new heavy compute; the point is to weigh competing axes honestly.

---

## TL;DR — the recommendation is regime-dependent

There is **no single best prime**; the answer splits by how deeply you concatenate and by
how much you trust an *un-modeled* physical-control cost.

| Regime | Recommended *p* | Representative code | γ | C | n |
|---|---|---|---|---|---|
| **Few-round, minimum single-round overhead** | moderate, **p ≈ 29–47** | `[[29,12,3]]₄₁`, `[[21,8,3]]₂₉` | 0.80, 0.88 | **3.2** | 29, 21 |
| **Balanced / realistic "sweet spot"** | **p ≈ 41–71** | `[[36,11,6]]₄₇`, `[[53,18,7]]₇₁` | 0.66, 0.56 | 4.7, 5.0 | 36, 53 |
| **Deep concatenation, control assumed cheap** | **largest constructible, p ≈ 89–97** | `[[67,22,9]]₈₉`, `[[73,24,9]]₉₇` | **0.51** | 6 | 67, 73 |
| **Below the RS threshold (p < 23)** | only if small-*p* control is *dramatically* cheaper | `[[234,55,5]]₁₇` | 0.90 | **38.9** | 234 |

**Headline tradeoff (p=23 → p=97, the full constructible RS ladder):** the overhead
*exponent* γ nearly halves (0.924 → 0.506) while the single-round cost C rises only
~1.5× (4.3 → 6.3) and the block size n rises ~4× (18 → 73). On **code-level metrics
alone, higher p is strictly better**. The thing that bounds *p* in practice is the
**un-modeled physical-control cost** (gate calibrations ~p², leakage channels ~p, state
prep), together with the **distance-computability wall** that freezes the constructible
family at the m=1 Reed–Solomon (RS) ladder. The biggest γ gains are in the **p ≈ 37–71**
band; beyond p≈71 the γ curve flattens (0.555→0.506 over Δp=26) while control cost keeps
growing ~p² — so the *marginal* return collapses and the practical sweet spot is
**p ≈ 41–71**, not the largest prime.

---

## 1. The cost framework (stated precisely)

Triorthogonal CSS code `[[n,k,d]]_p` consumes n noisy magic states at input error δ_in,
post-selects on the X-syndrome, and outputs k better states. Four figures of merit:

**(a) Overhead exponent γ** (deep-concatenation metric):
```
γ = log(n/k) / log(d)
```
γ<1 ("sublogarithmic") means the resource overhead to hit a target error grows *slower*
than the log of the target — the whole point of the search. Asymptotically the best
achievable yield is `γ₀(p)`, and it **falls with p like 1/ln p** (`qmsd.asymptotics.optimal_gamma`):

| p | 2 | 3 | 5 | 7 | 11 | 13 | 17 | 23 | 37 | 59 | 97 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **γ₀(p)** | 0.678 | 0.632 | 0.559 | 0.508 | 0.441 | 0.418 | 0.383 | 0.347 | 0.297 | 0.258 | 0.223 |
| **1/ln p** | 1.44 | 0.91 | 0.62 | 0.51 | 0.42 | 0.39 | 0.35 | 0.32 | 0.28 | 0.25 | 0.22 |

γ₀(p) tracks 1/ln p closely for p≳7. **More headroom below γ=1 the higher you go.** (The
optimal puncture fraction θ/(p−1) drifts slowly from 0.27 to 0.31 across this range.)

**(b) Single-round cost C** (`qmsd.distillation.cost / nbar_T`), at fixed δ_in:
```
n̄_T = (1 − (p−1)/p · δ_in)^n · k ,   C = n / n̄_T      (inputs consumed per accepted output)
```
C depends only on (n,k,p,δ_in); it is the *concrete* per-round price tag.

**(c) Per-round suppression** (`qmsd.distillation.delta_out_*`):
```
δ_out ≈ A_d / ((p−1) p^{d−1}) · δ_in^d        (per output; eq 38)
```
d sets the *order* of suppression; `A_d` (count of weight-d logical-Z operators) is the
prefactor. Lower A_d at fixed (n,k,d) = strictly better code.

**(d) Feasibility / constructibility — two gates:**
- **Singleton bound** for γ<1: need `D > 2√(p^m) − 2`. The first feasible m shrinks as p
  grows; at m=1 any p≥23 already clears it, while p<23 needs m≥2 (qutrit needs d≥7,
  effectively m≥7 — out of reach).
- **Distance-computability wall** (the load-bearing computational fact of this project):
  - **m=1 ⇒ the quantum code is *punctured MDS* ⇒ distance is CLOSED FORM**
    `d = r_max − k + 2`, `r_max = ⌊(p−2)/3⌋`, **independent of which columns are punctured.**
    No search. The entire p≥23 family is exactly enumerable.
  - **m≥2 ⇒ distance needs a meet-in-the-middle search** whose weight-3 enumeration over
    n=p² columns scales ~1e7 (p=23), 6.6e10 (p=97) — **infeasible past ~p=37.** So
    low-γ large-p codes at m≥2 *exist* but are **distance-uncomputable** with current tools.

The practical consequence: the **constructible frontier is exactly the m=1 RS ladder for
p≥23**, plus a few hand-certified m=2 points at p=17–19.

---

## 2. Per-axis analysis (with the actual numbers)

### Axis 1 — Overhead exponent γ: favors HIGH p

Min-γ RS code per prime (closed form, cross-checked in qmsd; `rs_family.py`):

| p | min-γ code | γ | p | min-γ code | γ | p | min-γ code | γ |
|---|---|---|---|---|---|---|---|---|
| 23 | `[[18,5,4]]` | 0.924 | 47 | `[[36,11,6]]` | 0.662 | 73 | `[[56,17,8]]` | 0.573 |
| 29 | `[[22,7,4]]` | 0.826 | 53 | `[[40,13,6]]` | 0.627 | 79 | `[[60,19,8]]` | 0.553 |
| 31 | `[[25,6,5]]` | 0.887 | 59 | `[[45,14,7]]` | 0.600 | 83 | `[[62,21,8]]` | 0.521 |
| 37 | `[[29,8,5]]` | 0.800 | 61 | `[[47,14,7]]` | 0.622 | 89 | `[[67,22,9]]` | 0.507 |
| 41 | `[[31,10,5]]` | 0.703 | 67 | `[[51,16,7]]` | 0.596 | 97 | `[[73,24,9]]` | 0.506 |
| 43 | `[[34,9,6]]` | 0.742 | 71 | `[[53,18,7]]` | 0.555 | | | |

**Diminishing returns are sharp.** Per-prime γ drops:
23→37: Δ0.124, 37→59: Δ0.200, 59→97: Δ0.094 — the steepest descent is the **p≈37–71
band**; past p≈71 the curve flattens (0.555→0.506 over Δp=26). So most of the reachable
γ-headroom is captured by p≈70, consistent with the 1/ln p shape (each *doubling* of p
buys a fixed Δγ≈ln2/(ln p)², which itself shrinks).

What γ buys: deep-concatenation resource cost scales ~ (log 1/ε)^γ. Halving γ from 0.92
to 0.51 at, e.g., log(1/ε)=30 is a `30^{0.92}/30^{0.51} ≈ 4×` resource reduction — and it
*grows* for deeper targets. **For deep concatenation the γ gain is real and dominant.**

### Axis 2 — Single-round cost C and block size n: grow with p (modestly)

From the per-prime (γ,C) Pareto frontiers (`rs_family.py` + Pareto enumeration):

| p | min-γ code | γ | **C** | **n** | cheapest γ<1 code (min C) | γ | **C** |
|---|---|---|---|---|---|---|---|
| 23 | `[[18,5,4]]` | 0.924 | 4.3 | 18 | `[[17,6,3]]` | 0.948 | **3.3** |
| 29 | `[[22,7,4]]` | 0.826 | 3.9 | 22 | `[[21,8,3]]` | 0.878 | **3.2** |
| 41 | `[[31,10,5]]` | 0.703 | 4.2 | 31 | `[[29,12,3]]` | 0.803 | **3.2** |
| 47 | `[[36,11,6]]` | 0.662 | 4.7 | 36 | `[[33,14,3]]` | 0.780 | 3.3 |
| 71 | `[[53,18,7]]` | 0.555 | 5.0 | 53 | `[[49,22,3]]` | 0.729 | 3.6 |
| 97 | `[[73,24,9]]` | 0.506 | 6.3 | 73 | `[[67,30,3]]` | 0.731 | 4.4 |

- **C grows slowly:** min-γ C goes 4.3 (p=23) → 6.3 (p=97), only ~1.5× across a 4× span of p.
- **The absolute-cheapest single round is ~C≈3.2**, hit by the **d=3, high-k codes around
  p=29–47** (`[[21,8,3]]₂₉`, `[[29,12,3]]₄₁`). The cheap end is *flat in p* — you don't
  need a large prime to get a cheap round.
- **n grows ~linearly with p** (18→73, ~4×). For m=1, n = p − k ≲ p, so block size is
  bounded by the prime itself — these are *tiny* codes (n<100 even at p=97). This is the RS
  family's quiet superpower: low-γ AND small blocks, simultaneously.

**Trade for the few-round regime:** you cannot amortize γ over many levels, so C and n
dominate. The (γ,C) frontier knee sits at **p≈41–71** (γ 0.55–0.70 at C 4–5); the rock-
bottom-C corner sits at **p≈29–47** (C≈3.2, γ 0.78–0.88).

### Axis 3 — A_d / suppression-per-round: the subtle axis (read carefully)

For the MDS quantum codes the prefactor is a **verified closed form**
`A_d = C(n,d)·(p−1)` (checked against `qmsd.weightdist` for p=23,29,31,37 — exact match:
67320, 204820, 1593900, 4275180). The **raw A_d explodes** with p:

| p | code | d | **A_d** | δ_out @ δ_in=0.01 |
|---|---|---|---|---|
| 23 | `[[18,5,4]]` | 4 | 6.7×10⁴ | 5.2×10⁻¹⁰ |
| 37 | `[[29,8,5]]` | 5 | 4.3×10⁶ | 8.3×10⁻¹³ |
| 47 | `[[36,11,6]]` | 6 | 9.0×10⁷ | 8.2×10⁻¹⁶ |
| 59 | `[[45,14,7]]` | 7 | 2.6×10⁹ | 8.2×10⁻¹⁹ |
| 73 | `[[56,17,8]]` | 8 | 1.0×10¹¹ | 8.2×10⁻²² |
| 97 | `[[73,24,9]]` | 9 | 9.3×10¹² | 5.6×10⁻²⁵ |

**The crucial honesty correction.** Naively "A_d is large ⇒ large p is bad." That is
*wrong as stated*. The raw A_d grows ~10⁸× from p=23 to p=97, but the **actual per-round
output error δ_out *improves* by ~15 orders of magnitude** — because the larger distance d
and the p^{d−1} denominator overwhelm the A_d prefactor. The effective coefficient
`A_d/((p−1)p^{d−1})` actually *falls* (0.25 at p=23 → 1.2×10⁻⁵ at p=97). **On absolute
per-round suppression, higher p wins decisively, not loses.**

So in what sense is MDS-A_d a genuine downside? Two precise senses, both narrower than
"large p is bad":

1. **Lost optimization freedom.** MDS distance/A_d are *puncture-invariant* — you are
   *stuck* with the combinatorial-maximum A_d = C(n,d)(p−1) and have **zero lever** to
   lower it. Contrast the qutrit cap-set result: at *fixed* `[[206,37,4]]₃` a cleverer
   puncture set cut A_d 880→**572** (1.54× per round, compounding to **8.6×** at 2 rounds,
   **8500×** at 3 rounds — `QUTRIT_PARETO.md`, `STRUCTURED_AD.md`). That lever **does not
   exist** for MDS codes.
2. **Where the lost lever matters: small p / small d / many rounds.** A forfeited constant
   factor compounds over concatenation levels. At large p the per-round suppression is
   already so steep (d=8–9) that one round overshoots almost any target and you run *few*
   rounds — so the lost constant is immaterial. At small p (d=3–4) you run *many* rounds
   and the constant compounds — exactly the qutrit regime where A_d-optimization is the
   *only* reachable win (since γ<1 is out of reach there at all; `QUTRIT_PARETO.md §1`).

**Net:** Axis 3 does **not** push toward small p for suppression quality. It says: large-p
MDS gives excellent but *un-tunable* suppression; small-p non-MDS gives weaker but
*optimizable* suppression, and that optimization is the only game in town when γ<1 is
unreachable.

### Axis 4 — Distance-computability / constructibility: a HUGE advantage for p≥23

This is decisive and easy to undervalue:

- **p≥23, m=1:** distance is closed form (MDS), A_d is closed form `C(n,d)(p−1)`. The
  *entire* γ<1 family — per-prime 2–7 Pareto-optimal codes, **21 codes** on the cross-prime
  (γ,C) envelope and **26** on the (γ,n) envelope — is enumerable in milliseconds with
  **no search and no distance bottleneck**. You can *certify* every code you propose.
- **p<23 needs m≥2**, where distance requires the MITM search that explodes past ~p=37.
  The one hand-found example, `[[234,55,5]]₁₇` (γ=0.90, the first certified γ<1 search code
  below the RS threshold), required a dedicated m=2 search and has **C=38.9** and **n=234**
  — i.e. ~9× the C and ~13× the n of `[[18,5,4]]₂₃`. It is, on every code-level axis,
  **dominated by `[[22,7,4]]₂₉`** (γ 0.826<0.90, C 3.9<38.9, n 22<234). It is a research
  milestone, not an operating point.
- **Large-p, m≥2 codes** would have even *lower* γ than the RS ladder, but they are
  **distance-uncomputable** today (past p≈37). So the constructible/certifiable frontier is
  *capped* at the m=1 RS ladder regardless of what exists in principle.

Constructibility therefore **reinforces** the p≥23 RS regime and **caps** the useful upper
end at whatever m=1 reaches (p≈97 here; the closed form runs to arbitrary p, but see Axis 5).

### Axis 5 — Physical control difficulty: the missing cost that bounds *p* (qualitative)

Everything above is a *code-level* accounting. It omits the dominant real-world cost: a
dimension-p qudit is physically harder than a qubit, and **monotonically harder in p**:

- **Two-qudit entangling gates** need ~O(p²) independent calibrations / pulse parameters
  (the Clifford generators act on a p²-dim two-qudit space); native high-fidelity p-ary
  gates are far less mature than qubit gates.
- **Leakage / error channels** multiply: p−1 excited levels means more leakage paths, more
  distinct depolarizing components, and a worse native δ_in than a comparable qubit — and
  the cost model takes δ_in as a *fixed input*, hiding this.
- **State preparation & readout** of p-level states is lower-fidelity and slower, which
  directly *raises* the effective δ_in that feeds C and δ_out.
- **The magic state itself** (a p-level non-stabilizer state) is harder to inject at higher p.

None of this is in γ, C, A_d, or δ_out. Crucially, this cost is **plausibly super-linear
in p** (gate-calibration ~p², channel count ~p), while the *marginal γ benefit is
sub-linear and decaying* (Δγ per Δp ∝ 1/(p ln²p), Axis 1). **Two opposing curves — a
decaying code benefit and a growing physical cost — generically produce an interior
optimum**, i.e. a finite sweet-spot p, not "as large as possible." We cannot locate it
numerically because the physical-cost curve is un-modeled; but its existence and the fact
that it sits at *finite* p is robust to the model's details.

---

## 3. Where the cross-prime Pareto envelope says to operate

Enumerating all γ<1 RS codes over p∈{23,…,97} and taking the non-dominated set:

**(γ, C) envelope — 21 codes.** From cheapest to lowest-γ:
`[[29,12,3]]₄₁` (γ0.80,C3.2) → `[[33,14,3]]₄₇` (0.78,3.3) → `[[34,13,4]]₄₇` (0.69,3.7) →
`[[38,15,4]]₅₃` (0.67,3.7) → `[[42,17,4]]₅₉` (0.65,3.7) → `[[51,20,5]]₇₁` (0.58,4.2) →
`[[59,24,5]]₈₃` (0.56,4.4) → `[[63,26,5]]₈₉` (0.55,4.5) → … → `[[73,24,9]]₉₇` (0.51,6.3).

**(γ, n) envelope — 26 codes**, same shape, n from 17 to 73.

Reading the envelope by regime:

**(a) Deep-concatenation regime (γ dominates).** Operate at the **low-γ end: p≈83–97**,
`[[67,22,9]]₈₉` / `[[73,24,9]]₉₇`, γ≈0.506, C≈6, n≈70. You pay ~1.5× the single-round C of
the small-p codes and get a ~2× lower overhead *exponent*, which dominates once you
concatenate. **But:** past p≈71 you are buying almost no further γ (0.555→0.506) for steadily
more physical control — so even in this regime the *practical* recommendation is the
**knee at p≈71–89**, not the literal maximum.

**(b) Few-round, low-overhead regime (C and n dominate).** Operate at the **cheap end:
p≈29–47**, the d=3–4 high-k codes (`[[21,8,3]]₂₉`, `[[29,12,3]]₄₁`, `[[34,13,4]]₄₇`),
C≈3.2–3.7. You give up γ (0.69–0.88) but you only run 1–2 rounds, so the asymptotic
exponent is irrelevant and the per-round price is what you pay. Going to p=97 here would
*raise* C by ~2× for a γ you never cash in.

**(c) The balanced sweet spot.** If you want most of the γ headroom *and* small, cheap,
certifiable codes *and* you respect that physical control grows with p, operate at
**p≈41–71**: `[[36,11,6]]₄₇` (γ0.66,C4.7,n36), `[[40,13,6]]₅₃` (0.63,4.6,40),
`[[53,18,7]]₇₁` (0.56,5.0,53). This band captures **~80% of the γ-drop from p=23 to p=97**
(0.924→~0.56 of the available 0.42) at **half the block size** of the p=97 codes and on the
flat-cheap part of the C curve — i.e. *before* both the γ curve flattens and the (assumed)
control cost runs away.

---

## 4. Recommendation

1. **Use a prime p ≥ 23 with m=1 Reed–Solomon codes.** This is the only regime that is
   simultaneously (i) γ<1, (ii) small-block (n≲p), (iii) cheap (C≈3–6), and (iv) *exactly
   constructible/certifiable with no search*. p<23 forfeits (i)+(iv) and, where γ<1 is even
   reachable (p=17 m=2), pays ~9× in C and is dominated by p=29.

2. **Pick the prime by concatenation depth:**
   - *Few rounds / modest target:* **p ≈ 29–47**, d=3–4, C≈3.2–3.7.
   - *Balanced default:* **p ≈ 41–71**, d=5–7, C≈4–5, γ≈0.56–0.70 — the recommended general-
     purpose choice.
   - *Deep concatenation, control cheap:* **p ≈ 83–97**, d=8–9, γ≈0.506, C≈6.

3. **Do not chase the largest prime.** The γ benefit decays like 1/ln p (and its
   *derivative* like 1/(p ln²p)); the per-prime γ gain past p≈71 is ≤0.05 while every
   physical-control cost keeps rising. The decaying-benefit/rising-cost crossing implies an
   **interior optimum at finite p** — our code-level data plus the qualitative control curve
   put it around **p≈40–70**.

4. **The qutrit/small-p lane is a *suppression-quality* play, not an overhead play.** If a
   platform's control is so much easier at p=3,5,7 that it dominates everything, you cannot
   get γ<1 there — but you *can* still win via **A_d-optimization** at fixed parameters
   (the verified `[[206,37,4]]₃` A_d 880→572, compounding 8500× over 3 rounds). That is the
   honest small-p recommendation: optimize the prefactor you *can* move, since the exponent
   is locked.

---

## 5. What would change this answer

The recommendation rests on assumptions that, if modeled, could move the sweet spot:

- **The missing physical-cost model is the whole ballgame for the upper bound.** We assert
  control cost grows ~p² (gates) to ~p (channels) and that it pushes the optimum to finite
  p, but we **do not have a quantitative model**. A platform where p-ary gates are *as cheap
  as* qubit gates (e.g. native spin-S or photonic qudits with cheap high-d operations) would
  push the optimum all the way to the largest constructible p (≈97 here, or higher via the
  m=1 closed form). A platform where cost explodes past p≈10 would invalidate the entire
  p≥23 regime and force the qutrit A_d-optimization lane. **Every "high p wins" statement
  above is conditional on this un-modeled cost being sub-dominant to the γ gain.**

- **δ_in is taken as a fixed, p-independent input.** Realistically δ_in *worsens* with p
  (more leakage, harder state prep). Since C and δ_out are steep in δ_in, a p-dependent
  δ_in(p) would penalize large p and pull the sweet spot down — a first-order correction the
  current model omits.

- **The distance-computability wall caps the *constructible* frontier, not the *existent*
  one.** Low-γ codes at large p, m≥2 exist (and would beat the RS ladder on γ) but are
  uncertifiable past p≈37 with current MITM tooling. A better distance/A_d algorithm
  (e.g. the structured-flat enumerator generalized off qutrit, or a distributed weightdist
  engine — see `STRUCTURED_AD.md`, `DISTRIBUTED_WEIGHTDIST_DESIGN.md`) would *open new
  operating points below the RS γ-floor* and could shift the deep-concatenation
  recommendation. Until then, "operate on the m=1 RS ladder" is a statement about *what we
  can build and certify*, not about what is optimal in principle.

- **γ is an asymptotic exponent; the few-round recommendation uses a single-round C.**
  A genuine multi-round total-overhead optimizer (compounding C and δ_out across a concrete
  level schedule to a fixed target ε) would sharpen the (a)/(b) boundary and might select
  *mixed-p* concatenations (cheap high-k codes early, low-γ codes deep). We have the
  per-round pieces (`qmsd.distillation`) but did not run the full schedule optimization;
  doing so is the natural next quantitative step.

- **A_d closed form is verified only for the MDS (m=1) family** (p=23,29,31,37, exact). The
  qualitative "A_d is un-tunable for MDS" claim is structural (puncture-invariance) and
  solid; the quantitative δ_out advantages of large p assume the eq-38/39 averaged model in
  `qmsd/distillation.py` is the right operational metric.
