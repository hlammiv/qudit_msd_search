# Beyond γ: Clifford cost, block error, and detect-and-correct rounds

*2026-11-01. Companion to `LEVEL_SETS.md` (cost section) and `PERFORMANCE.md`. Scripts in
`scratchpad/`, data in `runs/coverage/` and `runs/rm_m4/`.*

## 1. Clifford cost per round (`scratchpad/clifford_cost.py`, `runs/coverage/clifford_cost.json`)

Gate model, an upper bound with no circuit optimisation: after the transversal T layer, each X-type
check (row of G0) costs one two-qudit gate per nonzero entry, and each output is decoded with one gate
per nonzero of its logical row. G0 is row-reduced first and the logical rows are reduced against its
pivots; a minimum-weight basis would be cheaper but is NP-hard to find. Codes too large to row-reduce
here use the dense-fraction law nnz ≈ r(n − r)(1 − 1/p), which the row-reduced codes follow to within
15 %.

| code | checks | mean check weight | two-qudit gates per round | per output | method |
|---|---|---|---|---|---|
| [[206,37,4]]_3 | 14 | 130 | 6,413 | **173** | row-reduced |
| [[519,106,5]]_5 | 16 | 401 | 48,713 | **460** | row-reduced |
| [[7,1,3]]_7 cyc-B | 2 | 6 | 17 | **17** | row-reduced |
| [[1745,271,7]]_7 | 55 | 1206 | 455,805 | **1,682** | row-reduced |
| [[1746,270,6]]_7 | 56 | 1227 | 456,412 | **1,690** | row-reduced |
| [[1756,260,7]]_7 | 66 | 1268 | 457,274 | **1,759** | row-reduced |
| [[2093,308,5]]_7 | 18 | 1782 | 579,863 | **1,883** | row-reduced |
| [[9,2,3]]_11 | 2 | 8 | 30 | **15** | row-reduced |
| [[100,21,4]]_11 | 7 | 86 | 2,301 | **110** | row-reduced |
| [[13332,1309,22]]_11 | 1011 | 6278 | 20,997,783 | **16,041** | row-reduced |
| [[13,4,3]]_17 | 2 | 12 | 68 | **17** | row-reduced |
| [[237,52,6]]_17 | 14 | 210 | 13,500 | **260** | row-reduced |
| [[78455,5066,≥68]]_17 | 7304 | 66967 | 828,372,358 | **163,516** | dense-fraction estimate |
| [[15,4,3]]_19 | 2 | 14 | 80 | **20** | row-reduced |
| [[296,65,6]]_19 | 13 | 269 | 20,513 | **316** | row-reduced |
| [[123284,7037,≥57]]_19 | 10233 | 107102 | 1,849,643,593 | **262,846** | dense-fraction estimate |
| [[22,7,4]]_29 | 3 | 20 | 193 | **28** | row-reduced |
| [[666,175,6]]_29 | 15 | 628 | 117,789 | **673** | row-reduced |

The large codes' checks are essentially dense: a row-reduced check touches about two thirds of the
block. The magic-state saving of the large p = 7 codes (4.1–4.3× proved) therefore holds only while a
raw magic state costs more than a few hundred two-qudit Clifford gates; with cultivated or injected
inputs that is plausible, but it is an assumption every cost claim in `PERFORMANCE.md` makes and
should state. A sparse basis for the checks (lower-weight generators of the same row space) is the
obvious open question: the level-set and weighted codes have structured duals, so short generators may
exist, and they would cut this cost by the same factor.

## 2. Block error versus per-output error (`runs/weighted/block_error_p7.json`)

From the complete weight-d censuses of the four large p = 7 codes; leading order in ε, block error =
A_d (ε/(p−1))^d, per-output error = L_d (ε/(p−1))^d.

| code | weight-d accepted vectors A_d | per-output count L_d | outputs wrong per failing block | block ÷ per-output |
|---|---|---|---|---|
| [[1746,270,7]]₇ (headline, seed 19) | 276 | 193 | 189 of 270 | 1.4 |
| [[1756,260,7]]₇ | 24 | 1.29 | 14 of 260 | 18.6 |
| [[1746,270,6]]₇ | 24 | 1.33 | 15 of 270 | 18.0 |
| [[2093,308,5]]₇ | 1560 | 1338 | 264 of 308 | 1.2 |

Failures come in clusters, from 14 outputs to most of the block. Per-output error is the right input
to the next round only if each next-round block draws its inputs from different blocks of this round;
if a factory feeds one block's outputs into one next-round block, the block error governs, and it is
up to 19× worse. Aanha's Task B (`docs/aanha_mc_guide2`) measures this end to end.

## 3. Detect-and-correct rounds (`qmsd/rigorous_curve.py: CorrectingBoundCurve`)

A Z power is a Pauli, so instead of discarding every block with a nonzero syndrome one can decode to
the lowest-weight error of weight ≤ t consistent with the syndrome (unique while 2t < d), undo it and
keep the block. Rigorous, census-free consequences:
- acceptance ≥ P(at most t input errors), against (1 − ε)ⁿ for detection;
- a kept block fails only if the true error is within distance t of a nonzero logical kernel word, so
  the suppression order drops from d to d − t. Two bounds, the smaller is used: P(wt(e) ≥ d − t), and
  Σ_{w≥d} C(n,w) p^{w−d+1} Σ_{s≤t} C(w,s)(ε/(p−1))^{w−s} (Singleton count of kernel words).

**Monte Carlo of the t = 1 lookup decoder** (`scratchpad/mc_correct_t1.py`; the syndrome must be a
multiple of one column of G0, 10 470 keys for the p = 7 code, none ambiguous):

| code, ε | acceptance, detect → correct | logical failures (kept blocks) |
|---|---|---|
| [[1745,271,7]]₇, 3·10⁻³ | 0.51 % → 3.25 % (theory 0.53 % → 3.30 %) | 0 in 6 502 |
| [[100,21,4]]₁₁, 5·10⁻³ | 60.6 % → 91.0 % | detect 0; correct: 8 outputs in 2 blocks |
| [[100,21,4]]₁₁, 10⁻² | 36.6 % → 73.6 % | detect 0; correct: per-output 1.3·10⁻⁵ |
| [[100,21,4]]₁₁, 2·10⁻² | 13.3 % → 40.3 % | detect 0; correct: per-output 6.2·10⁻⁵ |

The decoder raises acceptance exactly as predicted and costs one order of suppression: the distance-4
code's failures under correction grow as about ε⁴ across the range (the proof guarantees order ≥ 3),
while pure detection shows none at these shot counts.

**Protocol costs with correction allowed** (`scratchpad/correcting_cost.py`,
`runs/rm_m4/correcting_costs.json`; every t from 0 to ⌊(d−1)/2⌋ offered to the optimizer, baseline =
standard library + Reed–Solomon family, rigorous bounds throughout):

| code | 1e-2 → 1e-30 detect / correct | 1e-3 → 1e-30 detect / correct | chosen t |
|---|---|---|---|
| [[13332,1309,22]]₁₁ | 1.16× / 1.17× | 1.13× / 1.18× | 6–8 |
| [[78455,5066,≥68]]₁₇ | unused / unused | unused / unused | — |
| [[123284,7037,≥57]]₁₉ | unused / unused | unused / unused | — |

**Correction does not rescue the large codes' cost.** It widens the window in which they accept, but
at the input errors where they are used (10⁻⁵ to 10⁻⁴) acceptance is already near 1, and reaching that
window is what costs; the cheap Reed–Solomon codes at p ≥ 17 still finish first. The idea is sound and
cheap (t = 1 is a table lookup), it is worth keeping for codes with distance to spare, and it is the
natural ingredient of a graded-output protocol (label kept blocks by correction weight, re-distil the
heavily corrected ones), which is not yet modelled.

## 4. Graded rounds with recycling (`qmsd/graded.py`, `scratchpad/graded_cost.py`, `runs/coverage/graded_cost_p7.json`)

A graded round decodes each block to its lowest-weight explanation (weight ≤ t), labels the block with
that weight (its *grade*), sends grades ≤ g forward, recycles grades g+1..t as inputs to the same round,
and discards the rest. For the small p = 7 library codes the grade statistics are **exact**: a dynamic
program over the qudits gives the joint distribution of syndrome and logical error on F_p^(r+k) (up to
8·10⁵ states for [[28,1,5]]₇); grade 0 reproduces the exact engine to every printed digit. Recycling is
a steady state in which recycled outputs join the inputs at their average error rate (offered for k = 1
codes only, where this is accurate). The optimizer sees every (t, g, recycle) choice as a separate code.

| ε₀ → target | library | library + graded | library + large codes | … + graded |
|---|---|---|---|---|
| 1e-2 → 1e-10 | 37.1 | 36.7 (1.01×, picks [[28,1,5]]₇ graded) | 37.1 | 36.7 |
| 1e-2 → 1e-15 | 52.6 | 52.6 | 48.4 | 48.4 |
| 1e-2 → 1e-30 | 210.3 | 210.3 | 48.4 | 48.4 |
| 3e-3 → 1e-15 | 31.6 | 31.6 | 31.6 | 31.6 |
| 1e-3 → 1e-15 | 28.8 | 28.8 | 28.8 | 28.8 |
| 1e-3 → 1e-30 | 141.0 | 141.0 | 39.0 | 39.0 |

**Grading saves at most about 1 %, for two structural reasons.**
- *Distance-3 codes cannot recycle.* A corrected "single error" is too often a misread double error: at
  ε = 1 % the grade-1 blocks of [[7,1,3]]₇ carry 2.5 % error, worse than the raw inputs, so the optimizer
  declines. Grading needs d ≥ 4 and lower input error.
- *A recycled output replaces one of n inputs.* With k = 1 the saving is bounded by k·P(recycle)/n,
  about 1/28 × 0.25 ≈ 1 % for [[28,1,5]]₇ at ε = 1 %, which is what the optimizer finds. High-rate codes
  would recycle more per block, but their failures are clustered (section 2), so their corrected blocks
  carry many correlated output errors at once.

Conclusion: syndrome-dependent grading and recycling is cheap and harmless but not worth building
into a factory for these codes. Together with section 3, error correction inside distillation rounds is
closed as a cost lever here; the remaining levers are cheaper Clifford circuits (sparse check bases)
and better small codes at each prime.

## 5. Sparse check bases and light logical operators (`scratchpad/sparse_basis.py`, `runs/coverage/sparse_*.json`)

Same gate model as section 1. Checks: structured light vectors of the check space (p = 7 weighted codes:
(g−1)(g−6) times products of three linear forms; level-set codes: g times hyperplane indicators,
weight ≈ |H \ O|), selected greedily for independence and completed with row-reduced rows. Outputs:
each logical operator H1_i + rowspace(G0) is pushed down by greedy descent (repeatedly add the
multiple of a check that lowers the weight most), on sampled outputs.

| code | check gates: row-reduced → sparse | mean logical weight: row-reduced → descended | gates per output: before → after |
|---|---|---|---|
| [[1745,271,7]]₇ | 66 316 → 58 188 (1.14×) | 1435 → 1403 of n = 1745 | 1680 → 1623 (**1.04×**) |
| [[13332,1309,22]]₁₁ | 6.35·10⁶ → 4.54·10⁶ (1.40×; 471 structured rows) | 11 196 → 11 000 of n = 13 332 | 16 045 → 14 469 (**1.11×**) |

**Decoding the outputs dominates (85–90 % of the gates) and its operators stay nearly dense** (about
0.8 n after descent). The lightest check-space vectors are not the lightest Reed–Muller words: a check
must vanish on every output, so the 294-point (p = 7) or 968-point (p = 11) RM words do not qualify,
and the earlier "up to 4–12×" estimate was wrong. For these dense codes the Clifford cost per output
per round is about (1 − 1/p)·n regardless of basis. The only route to a real reduction would be an
encoder that uses the Reed–Muller structure (evaluation of a polynomial axis by axis, m·p gates per
qudit, the analogue of the fast RM encoders); whether that is compatible with outputs at arbitrary
points (for the weighted codes) or on a cubic surface (level sets) is open.

**Bottom line for cost claims.** The large codes' magic-state savings (4.1–4.3× at p = 7) come with
a Clifford bill of about 1 600 two-qudit gates per output per round against about 17 for [[7,1,3]]₇;
the p = 11 level-set code needs about 14 500. A cost claim in the paper must state the price ratio
of a magic state to a Clifford gate at which the large codes break even.

## 6. One-shot rounds with an axis-by-axis encoder (`scratchpad/fast_encoder.py`, `scratchpad/oneshot_check.py`, `runs/coverage/fast_encoder.json`)

Section 5 concluded that no choice of basis makes the checks or the output decoding cheap. Both can be
avoided. A round only needs the state Σ_{f ∈ RM_p(r,m)} |f(x)⟩ on the whole grid F_pᵐ: apply the noisy
T on every non-output qudit, measure those qudits in the X basis (the outcomes give the syndrome
classically and teleport the result), and the output qudits are left in U⁻¹|+⟩ᵏ up to a Pauli frame,
because Σ_x f(x)³ = 0 on the whole grid makes transversal U act trivially on the code state. There is no
check-measurement circuit and no per-output decoding.

That state is prepared by putting |+⟩ on the coefficient registers of the monomials of degree ≤ r,
|0⟩ on the rest, and applying the univariate evaluation map V[j,t] = jᵗ along every line of axis 1, then
axis 2, …, axis m (polynomial evaluation on a grid factorises as V⊗V⊗…⊗V). In stage i a line has only
the inputs whose exponent still fits under r, which prunes most of the work. Each line map is
implemented as an invertible p × p matrix and its two-qudit add-multiple gates are counted exactly by
Gauss–Jordan elimination (same convention as section 1: scalings and relabels free). Weighted supports
are prepared on the whole grid; points outside the support get no T (charge 0) and are measured with
the rest.

**Checks.** The factorisation equals direct evaluation and spans exactly rm_generator's space
(p, m, r) = (3,2,1), (5,2,2), (7,2,3), (5,3,3). The one-shot round was simulated exactly (state vectors
over all codewords): for RM(1,1) at p = 7 all 16 807 outcomes with nonzero probability leave the output in
Z^a U⁻¹|+⟩ with fidelity 1.000000000000; sampled outcomes (2000–3000 each) give the same for RM(3,1) at
p = 13, RM(1,2) at p = 7, and a two-output RM(1,2) code at p = 5.

| code | two-qudit gates per output per round: systematic encoder (section 1) | one-shot, axis-by-axis | ratio |
|---|---|---|---|
| [[206,37,4]]₃ | 173 | 29 | 6.0× |
| [[519,106,5]]₅ | 460 | 46 | 10.0× |
| [[1745,271,7]]₇ | 1 682 | 98 | **17×** |
| [[2093,308,5]]₇ | 1 883 | 86 | 22× |
| [[100,21,4]]₁₁ | 110 | 45 | 2.4× |
| [[13332,1309,22]]₁₁ | 16 041 | 233 | **69×** |
| [[237,52,6]]₁₇ | 260 | 73 | 3.6× |
| [[78455,5066,≥68]]₁₇ | 163 516 | 567 | **288×** |
| [[296,65,6]]₁₉ | 316 | 89 | 3.6× |
| [[123284,7037,≥57]]₁₉ | 262 846 | 701 | **375×** |
| [[666,175,6]]₂₉ | 673 | 114 | 5.9× |

Small codes ([[7,1,3]]₇ and the Reed–Solomon codes) cost 15–28 gates per output; the cyclic codes are not
grid-evaluation codes, so the method does not apply to them, and it does not need to.

**This reverses the conclusion of sections 1 and 5.** With the one-shot round the large codes cost
100–700 two-qudit gates per output per round, within a factor 4–40 of the small codes instead of
100–10 000, and per output of the *final* product the large codes are used once where small codes are
chained. The γ < 1 codes are no longer Clifford-prohibitive.

Caveats, all open:
- *Fault tolerance of the encoder.* A fault early in the transform spreads along a whole line (p
  qudits per stage, up to pᵐ⁻ⁱ after later stages). In distillation the Clifford circuit normally runs on
  error-corrected logical qudits, so such faults are suppressed by the outer code and the gate count is a
  space-time cost; if the circuit runs on physical qudits, its fault distance must be analysed
  separately (this is the "circuit-level distance" item, now the most important open metric).
- *Connectivity.* The stages act along lines of an m-dimensional grid; on planar hardware they need
  routing, which this count ignores (as did section 1).
- *All-to-all T.* The noisy T is applied to every non-output qudit, as before; inputs are teleported in.

## 7. Fault spreading in the one-shot encoder (`scratchpad/encoder_faults.py`, `runs/coverage/encoder_faults_*.json`)

Every single Pauli fault location of the axis-by-axis circuit (after each add-multiple gate, on both of its
qudits, plus preparation faults) was propagated exactly through the remaining linear gates, with the
circuit run on **physical** qudits. A final Z image z is *trivial* if it lies in C^⊥, *detected* if it
changes the X-measurement syndrome (G0 z_N ≠ 0), and *harmful* if it is undetected but changes the
outputs. An X image is equivalent, modulo a codeword shift, to a shift supported on the measured qudits;
after the T layer that acts like Z-type errors there, which are detected unless they form a kernel word
of weight ≥ d.

| code | fault locations | Z harmful (undetected) | of which: fault on an output after its last gate | from spreading | harmful locations per output |
|---|---|---|---|---|---|
| [[100,21,4]]₁₁ | 2 019 | 63 (3.1 %) | 51 | 12 | 3.0 |
| [[237,52,6]]₁₇ | 7 891 | 70 (0.9 %) | 68 | 2 | 1.3 |
| [[519,106,5]]₅ | 10 343 | 223 (2.2 %) | 182 | 41 | 2.1 |
| [[2093,308,5]]₇ | 55 519 | 811 (1.5 %) | 722 | 89 | 2.6 |

**Correction (2026-11-04): the table above counts Z faults only, and the X analysis below it was wrong.** An X
image was treated as equivalent, modulo a codeword shift, to a shift on the measured qudits and hence detected. That
fails: U acts on the inputs between the two, and an exact state-vector round shows an X on an output after its last
gate is always accepted and leaves fidelity 1/p (`tests/test_sdim_round.py`). Counting these too:

| code | Z harmful | X harmful | of which on an output directly | harmful locations per output (was) | input-escape term, Σ p^(−wt) (1 − 1/p) |
|---|---|---|---|---|---|
| [[100,21,4]]₁₁ | 63 | 118 | 117 | **8.6** (3.0) | 0.7 |
| [[237,52,6]]₁₇ | 70 | 448 | 448 | **10.0** (1.3) | 0.7 |
| [[519,106,5]]₅ | 223 | 361 | 346 | **5.5** (2.1) | 65 |
| [[2093,308,5]]₇ | 811 | 1 180 | 1 168 | **6.5** (2.6) | 36 |

The X count is higher because an output is the target of many gates, and an X on a target spreads only if that
qudit is later a control. The last column is what the model drops: an X image touching few inputs escapes with probability
about p^(−wt), since the random Z it becomes can be the identity. These location counts are for single-qudit Paulis of one
power. A two-qudit depolarizing fault is usually rejected by its other parts, so the floor per gate error is much
smaller; `qmsd/sdim_round.py` gets it jointly (about 0.004–0.12 p_gate on the codes tried).

**Spreading is not the problem.** Each output has 5–10 undetected first-order fault locations, almost all
of them faults on the output qudit itself, which every distillation circuit has; spreading
adds under 12 % on top (Z and X together). The encoder does not create the large correlated undetected errors one might fear,
because a fault early in the transform spreads to many qudits and is then caught by the checks.

**Acceptance is the cost.** The other 97–99 % of fault locations are detected and discard the round.
The circuit has about 26 locations per input qudit for [[2093,308,5]]₇ (55 519 against 2 093 inputs), so
the Clifford fault rate must be about 26× below the input magic-state error for the circuit not to dominate
the discard rate. That is the usual regime (Cliffords on error-corrected logical qudits, magic states noisy),
and the systematic encoder of section 1 would have about 17× more locations. On physical qudits with
comparable gate and input error rates the large codes would be throttled by acceptance, not by
undetected errors.

## 8. Combined cost: magic states and Clifford gates together (`scratchpad/total_cost.py`, `runs/coverage/total_cost.json`)

Objective per final output: R × (raw magic states) + (two-qudit Clifford gates), with R the cost of one magic
state in Clifford gates. Per round attempt a code costs n raw states and g gates (one-shot axis-by-axis count
for grid RM/RS codes, systematic encoder count for cyclic and doubled codes); per output of round i,
M_i = a_i M_{i−1} and Y_i = a_i Y_{i−1} + g_i/(k_i P_i) with a_i = n_i/(k_i P_i). The optimizer searches chains
of up to five rounds. Baselines: the standard library plus the whole Reed–Solomon family (exact curves). Large
codes: rigorous curves (p = 7: complete weight-7 censuses plus the weight-8 bound; p ≥ 11: census-free Singleton bounds).

Saving = cost(small codes only) / cost(with the large codes):

| p, ε₀ → target | R = 1 | R = 10 | R = 100 | R = 10⁴ | protocol with large codes |
|---|---|---|---|---|---|
| 7, 1e-3 → 1e-15 | 1.00 | 1.00 | 1.00 | 1.00 | small codes only |
| 7, 1e-3 → 1e-30 | **2.62×** | **3.23×** | **3.55×** | **3.60×** | [[6,1,?]]₇ RS → [[1746,270,7]]₇ (39.2 states, 171 gates per output) |
| 7, 1e-2 → 1e-30 | **3.24×** | **3.92×** | **4.27×** | **4.32×** | [[7,1,3]]₇ → [[1746,270,7]]₇ (48.7 states, 217 gates) |
| 11, 1e-3 → 1e-30 | 1.00 | 1.00 | 1.06× | 1.13× | [[8,3,2]]₁₁ → [[13332,1309,22]]₁₁ for R ≥ 100 |
| 17, 1e-3 → 1e-30 | 1.00 | 1.00 | 1.00 | 1.02× | [[13,4,3]]₁₇ → [[237,52,6]]₁₇ for R ≥ 1000 |
| 11, 17 at 1e-15 | 1.00 | 1.00 | 1.00 | 1.00 | small codes only |

**The p = 7 γ < 1 code wins at every price ratio**, by 2.6× even when a magic state costs one Clifford gate:
with the one-shot encoder it adds only ~200 gates per output, while the small-code chains it replaces need
several more rounds. At p ≥ 11 the Reed–Solomon family is already so good that the large codes add at most
13 %, and only when magic states are expensive. At 1e-15 targets no large code is used anywhere. (The p = 7
surrogate gives 2.6–3.6× at 1e-3 → 1e-30. Since 2026-11-11 the proved curve carries the rigorous weight-8 bound of
`WEIGHTED_SUPPORTS.md`, and the proved 1e-3 row equals the surrogate one; before it the row was 2.16 / 2.73 / 3.04 / 3.09× via
[[7,1,3]]₇ at 45.6 states and 209 gates. `scratchpad/total_cost_w8.py`, `runs/coverage/total_cost_w8.json`.)
