# D_PROOF_MAP — Triorthogonal no-go over small prime fields (honest proof map)

**Status banner.** The *literal global* conjecture is **FALSE** (verified counterexample, L8).
The *scoped m=2 / minimal-ambient* no-go is **OPEN**: its proposed crux (L7) was **REFUTED**,
and the genuine obstruction has been re-identified but not yet proven. What survives is a
**clean reduction with a corrected, still-open crux**, plus a rigorously calibrated p=17 boundary.

---

## 1. The conjecture, the levers, the anchor

**Conjecture (target).** For `p ≤ 13`, no triorthogonal CSS magic code over `F_p` achieves `γ < 1`.

**Column/moment setup.** X-stabilizer `G0` is `r × n_c` over `F_p`, `r = dim G0`. Columns are points
`x_i ∈ F_p^r`. Triorthogonality `⇔ Σ x_i⊗x_i = 0` and `Σ x_i⊗x_i⊗x_i = 0` (mod p) — this is exactly
`qmsd.is_triorthogonal` (ground truth, established in L1). Z-distance `d = min #` of F_p-dependent
columns (`min_dependent_columns`). With `k` logicals, `γ = ln(n_c/k)/ln(d)`, so **`γ<1 ⇔ d > n_c/k`.**

**Two rigorous levers.**
- **(R1) Singleton ceiling:** `d ≤ r+1`; defect `Δ = (r+1) − d ≥ 0`.
- **(R2) Cap floor:** `d ≥ 3` forces the columns to be a projective cap; `n_c > (p^r−1)/(p−1)` forces `d ≤ 2`.

**Forced-grid anchor (prior, rigorous + verified).** For evaluation/Veronese configs, triorthogonality
forces `1_S ∈ RM(t,m)` with `2t<p` ⇒ `S` = full grid ⇒ only ordinary RM puncturing ⇒ flat-capped
below `γ<1` for `p≤13`. The general conjecture *drops* the evaluation-code hypothesis.

**Boundary mandate.** A correct theorem must give `γ≥1` for `p≤13` yet *allow* the realized
`p≥17` crossings ([[237,52,6]]₁₇ γ=0.847; [[293,68,5]]₁₉ γ=0.908). A "proof" that also kills p≥17 is wrong.

---

## 2. The proof chain, lemma by lemma (VERIFIED status)

| ID | Content | Status |
|----|---------|--------|
| **L1** | Dictionary: code ⟷ projective arc (`d = 1 + arc strength`) ⟷ vanishing deg-2,3 Veronese moments; `is_triorthogonal ⇔ M2=M3=0`. | **PROVED** (3× SOUND) |
| **L2** | R1+R2 in geometric form; records that the **defect Δ is a dead lever** (Δ=0 reachable at p≤13; crossings have Δ=12,9,7,6, non-monotone). | **PROVED** (3× SOUND) |
| **L3** | Hadamard/Veronese reformulation `is_triorthogonal ⇔ 1⊥C*2 ∧ 1⊥C*3`; **negative**: only 2 functionals, no bound on `dim C*2`, so Schur-square rigidity never fires (measured `dimC*2` ≫ window). | **PROVED** (3× SOUND); reformulation only, no no-go content |
| **L4** | **Boundary engine.** Moment curve is MDS (`d=r+1`, Δ=0) and triorthogonal iff `3(r-1)<p-1`. ⇒ `d_max = ⌊(p+1)/3⌋+1` = 3,3,5,5 for p=5,7,11,13; first `d=6` triorthogonal-MDS sub-cap at **p=17** (`3·4=12<16`). Reproduces the boundary, does **not** kill p≥17. | **PROVED** (3× SOUND) |
| **L5** | No-logical on the MDS branch: a defect-0 NRC cap admits `k=0` logicals (so the high-d MDS route to γ<1 is closed). | **CONDITIONAL** — conclusion holds, *written proof has a gap* |
| **L6** | Line word pins distance (m=2): a line hit `L` times yields an explicit weight-`(d_RM−L)` dual word in `ker(G0)` ⇒ `d ≤ d_lines = d_RM − max_ℓ|S∩ℓ|`. UPPER bound only. | **PROVED** (3× SOUND) |
| **L7** | **CRUX (as stated):** length forces line-spread; combined with L4+L6 gives `d ≤ n_c/k` for p≤13 m=2. | **REFUTED** |
| **L8** | Scope + boundary-faithfulness gate: literal global claim FALSE; defect crux FALSE; governing quantity is `k` not Δ; scope is m=2/minimal-ambient; p=5 m≥4 escape and p≥17 m=2 crossings must survive. | **PROVED** (3× SOUND) |

### L5 — CONDITIONAL (gap in proof, conclusion true)
The written "CORE STEP" derives `P=0` only for **k=1**: it asserts `pair(L,L): Σ_nonmagic B² + P² = 0`,
but triorthogonality sums over **all** magic columns, so for `k≥2` one gets only `Σ_j P_j² = 0`, a sum
of squares that does **not** force each `P_j=0` over `F_p` (isotropic forms). Reviewers also flagged that
the all-ones row `R_0` used in the subtraction is in `rowspan(G0)` for the *unpunctured* NRC but **not**
once a logical is punctured in (verified false in every tested case). **However the conclusion `k=0`
survives**: the correct argument is the Gram identity `M·Mᵀ = 0` for the `k×k` logical-by-magic block
(from `pair(L_a,L_b) − triple(R_0,L_a,L_b)`); full rank needs `rank(M)=k`, but a totally-isotropic
subspace has `rank ≤ ⌊k/2⌋ < k`, contradiction. Reviewers brute-forced p=5,7,11,13 and found **no**
full-rank isotropic block. So L5 is **salvageable** — needs the Gram/isotropic rewrite to become PROVED.

### L7 — REFUTED (this is the heart, and it broke)
Read literally — a universal spread bound over valid triorthogonal m=2 puncture sets — L7 is **false**
for **p=11 and p=13**. Two independent fatal problems, both verified against ground truth:
1. **The 3rd-moment hypothesis is INERT.** `G0 = shorten(RM_p(r,2), S)` is a subspace of a triorthogonal
   space restricted to its support, so `is_triorthogonal(G0)=True` for *every* full-rank `S` (confirmed
   35/35 random S). L7 thus collapses to a hypothesis-free `(k;t)`-arc inequality
   `m_t(2,p) ≤ p²/(d_RM−t+1)`.
2. **That arc inequality is false.** Explicit p=13, k=34 set: `full_rank=True`, `is_triorthogonal=True`,
   `max_line=4 < d_RM − n_c/k = 5.03` — violated by 2; same at p=11. So `m_4(2,13) ≥ 34 > 28.17`.

Moreover L7 has the causation **backwards**: spread does not force large line-hits; spread **collapses
the rank** (`dim G0=2` in the counterexample), so the *true* `d=1 ≪ d_lines=5`. The real obstruction is
**R1 (Singleton/dim-G0 ceiling) + L4 (2D full-span codeword cap)**, not the planar-arc/Segre mechanism.

---

## 3. What is RIGOROUSLY ESTABLISHED end-to-end

**Not the conjecture.** The chain does **not** prove the global no-go, and indeed L8 *disproves* the
literal global statement.

What is rigorously nailed down:
- **(A) The literal global conjecture is FALSE** — `[[519,106,5]]₅` (p=5, m=4, dim G0=16, d=5 > n_c/k=4.896,
  γ=0.987), re-verified with ground-truth tools. The defect-crux ("γ<1 forces small defect") is also false
  (this code has Δ=12). [L8, PROVED.]
- **(B) The correct invariant is `d > n_c/k`, governed by `k`, not the Singleton defect.** [L1+L2+L8.]
- **(C) The p=17 boundary is reproduced rigorously and parameter-free** by the single inequality
  `3(r-1) < p-1`: triorthogonal-MDS distance is `d_max=⌊(p+1)/3⌋+1` (≤5 for p≤13), and a `d=6` MDS
  sub-cap first exists at p=17 — without killing p≥17. [L4, PROVED.]
- **(D) Two clean, correct distance levers:** the Singleton/cap geometry (R1,R2) [L2] and the explicit
  line-supported dual word `d ≤ d_RM − max_line` for m=2 [L6].

**Honest characterization:** a **verified reduction**. The scoped (m=2 / minimal-ambient) no-go is reduced
to a single crisp incidence/rank statement (Section 4), with all surrounding machinery (L1–L4, L6, L8)
proved and the L5 branch closed modulo a known rewrite. The originally-proposed crux (L7) is *refuted*, so
the reduction now points at a **corrected** target, not L7's arc bound.

---

## 4. The single cleanest remaining gap (the corrected crux)

> **Corrected CRUX (rank/2D form, m=2, p ≤ 13).** Let `S ⊆ AG(2,p)` be a full-rank puncture set of
> `RM_p(r,2)` with `3r < 2(p-1)` (so `G0` is automatically triorthogonal), `|S| = k`, `n_c = p² − k`.
> Prove that the **true** Z-distance satisfies `d ≤ n_c/k` for all `p ≤ 13`. Equivalently: to achieve
> `true d > n_c/k` a triorthogonal m=2 cap must be **rank-rich** (`dim G0 ≳ d−1`, by R1) **and** carry no
> weight-`≤5` two-dimensional (full-span) codeword (L4); show the puncture density forced by `k>0`
> logicals at `p≤13` makes a weight-`≤5` 2D codeword unavoidable, capping `d ≤ 5` and `d ≤ n_c/k`,
> while the same argument first fails at `p=17` (degree-5 NRC sub-cap, `3·5=15<16`).

This is the honest replacement for L7. It must **not** be phrased via the planar arc/Segre bound (that is
the refuted route) and must use the **true** governing pair: Singleton rank ceiling + 2D-codeword cap.
The empirical signal is exactly right (best m=2 γ: p5≈2.39, p7≈1.79, p11≈1.13, p13≈1.09, p17≈0.85),
so the target is to *prove* the `d ≤ n_c/k` crossing-blocker for p≤13 m=2 — currently **OPEN**.

---

## 5. Is the conjecture strongly supported, and is the p≥17 boundary respected?

- **Literal global form:** refuted — do **not** support it.
- **Scoped m=2 / minimal-ambient form:** **strongly empirically supported.** Across exhaustive and
  adversarial searches at p∈{5,7,11,13}, m=2, every triorthogonal cap has `d≤5` and `d≤n_c/k` (zero
  crossings; best γ≈1.06–1.13 at p=13/11). The proved pieces (R1, L4 d-ceiling, L6 line word) all point
  the same way; only the *quantitative closure* (Section 4) is missing.
- **p≥17 boundary:** **respected and rigorously calibrated.** L4's `3(r-1)<p-1` lands the first `d=6`
  MDS sub-cap exactly at p=17; L8's gate confirms [[237,52,6]]₁₇ and [[293,68,5]]₁₉ are classified γ<1.
  No proved lemma over-kills p≥17 (every adversarial counterexample search confirmed this per-lemma).

---

## 6. Concrete next steps to finish

1. **Rewrite L5 to PROVED** via the Gram/isotropic argument: from `pair(L_a,L_b) = triple(R_0,L_a,L_b)`
   derive `M Mᵀ = 0` for the logical-by-magic block, then `rank(M) ≤ ⌊k/2⌋ < k` kills full rank. Quick win;
   reviewers already found no full-rank isotropic block at p=5,7,11,13.
2. **Attack the corrected crux (Section 4) as a rank/2D statement, not an arc bound.** Concretely:
   prove that at p≤13, any full-rank m=2 puncture set with `dim G0 ≥ 5` (needed for `d=6` via R1) contains
   a 2D full-span codeword of weight `≤5` in `ker(G0)`, forcing `d≤5≤n_c/k`. Use the GRM weight hierarchy
   (Delsarte–Goethals–MacWilliams α=1) for the *second*-lowest weight class, beyond L6's line words.
3. **Pin the exact scope** ("minimal-ambient / no-smaller-code") into a precise hypothesis; check m=3
   separately (prior data: plane/2-flat words cap d at 4 then 2, γ≥1.21 — likely easier).
4. **Boundary regression test:** any candidate crux lemma must be auto-checked against the L8 gate
   (p=5 m=4 escape intact; [[237,52,6]]₁₇, [[293,68,5]]₁₉ classified γ<1) before acceptance.
5. **Close the small loose end** `m_t(2,7)` exact arc numbers — only relevant if one still wants the
   (now-abandoned) arc framing for documentation; not on the critical path.
