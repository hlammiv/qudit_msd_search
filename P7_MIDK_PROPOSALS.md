# p=7 m=4 mid-k window — synthesized direction dossier

*Multi-agent workflow (4 proposal lenses: algorithms / nogo / existence / reformulate, each adversarially vetted, then synthesized). Drop-in companion to `D_P7M4_WINDOW.md`. All figures re-verified against the repo (2026-07-31).*

**Verified against the repo:** `d_RM(16,4)=21`; plane `RM_7(4,2)=[49,15,21]`; hyperplane `RM_7(10,3)=[343,226,21]`; 139650 2-flats / 2800 hyperplanes in AG(4,7); int64 overflow bites at `dim G0 ≥ 23` i.e. `k ≤ 303` (so `k≥304` is int64-safe); coplanarity cap `d≤18` closes exactly `k≤126` (fails at 127); `weightcount` hard-caps at `d≤6` (line 105); `mindist` raises on `p**r > int64` (line 46); `structured_pe` rejects `RM_7(4,2)` (beta=−2); `stern_isd.py` is F_3-hardcoded throughout; `scratchpad/p7m4_framework.py` does **not** exist.

## Bottom line up front

The window **cannot be rigorously CLOSED with any tool in this repo** without proving the *same* open crux the whole Direction-D program already faces at m=2 — the **full-span DGM second-weight-class cap** (`D_CRUX_REDUCTION.md` §4, unproven even for p=11 m=2). Every purely geometric / averaging / LP / GHW route **provably tops out above `needed(k)`** because the binding codewords are large full-span words swallowed by `S`, which flat structure cannot see (validated: at k=314 the min-weight bound gives 8, true d=3). What *is* achievable: (1) **tighten both edges rigorously** to shrink the certified-open interval from `[110,312]` to `[127,303]` (cheap, sure, machine-verifiable), and (2) run a **one-sided F_7 Stern sweep** that per-S-closes the high-k sub-band down to ~250 and hunts for a survivor in the mid band — evidence, never a proof.

## The map, corrected

| Region | Status | Mechanism |
|---|---|---|
| `k ∈ [110,126]` | **closable now (A1)** | coplanarity floor `d≤18 < needed` — extends repo's `[110,114]` by 12 |
| `k ∈ [127,~250]` | **genuinely open** | true `d ≈ 9–18`; no certifier reaches it; where a γ<1 survivor would hide; Stern infeasible at the required weight |
| `k ∈ [~250,303]` | **per-S closable (B1)** | collapse band, true `d ≈ 7–9`; F_7 Stern finds low-weight full-span witnesses |
| `k ∈ [304,312]` | **closable now (A2)** | int64-safe (dim G0 ≤ 22), `min_dependent_columns(d_max=4)` → d=2–3 |
| `k ≥ 313` | already rigorous | Hamming/MDS redundancy bound |

Net after A1+A2: certified-open interval collapses to **`[127,303]`**; the intrinsically-hard core (needs the crux) is **`[127,~250]`**.

---

## Part A — CLOSE IT (no-go direction — the repo's lean)

### A1 · Coplanarity floor + universal averaging cap — **RANK 1 (do this)**
*(merges: nogo-P1 "coplanarity floor", existence-P1 "invariant-LP cap", reformulate-P2 salvage — three lenses independently derived the same lemma)*

- **Idea.** Any `|S|≥3` contains 3 points that share a 2-flat; over F_7 a plane has 8 parallel line-classes, so some class puts the 3 points on 3 distinct parallel lines, which are the support of a genuine weight-21 min-weight codeword of `RM_7(16,4)`. That word loses ≥3 points to `S`, so `d(S) ≤ 18` **for every S**. Since `18 < needed(k)` exactly for `k ≤ 126`, this closes `[110,126]` (γ≥1). The invariant-LP cap `d ≤ ⌊21(2401−k)/2401⌋` is the same averaging argument but strictly weaker (closes only `[110,120]`) — keep it as the "averaging provably can't reach the bulk" corollary.
- **Blocker beaten.** All three, at the edge: it computes **no distance** — pure incidence geometry, immune to overflow / d≤6 cap / large-d wall.
- **Cost / where.** Trivial, **local, minutes, analytic.** Crossover is k=126 (n/k=18.056) → 127 (17.906).
- **First step.** Lift `structured_m3`'s plane-certificate emitter to a 2-flat of AG(4,7); emit one explicit weight-18 codeword for `k=124`, assert `X_stab·v ≡ 0 (mod 7)`. Fold the lemma into `D_P7M4_WINDOW.md`, replacing the loose "~k/7 per hyperplane" heuristic.
- **Honest ceiling.** 4 points need not be coplanar, so `d≤18` is the flat route's floor — it cannot touch `[127,303]`. A sure tightening, not a resolution.

### A2 · Run the int64-safe MITM sweep on `[304,312]` — **RANK 2 (do this)**
- **Idea.** `dim G0 = 326−k ≤ 22` for `k≥304`, so the encoder `7^dim = 3.9e18 < int64` does **not** overflow. `min_dependent_columns(d_max=4)` certifies the true (small) d there; the doc reports the collapse (k=314→3, k=320→2) but never ran a certified sweep of `[304,312]` itself.
- **Blocker beaten.** #1 doesn't fire in this range; #2 irrelevant (d is 2–3, well under 6).
- **Cost / where.** **Local, minutes-to-~1h total** (9 values of k). Run under `ulimit -v 4GB` — the k=305 build transiently needs ~0.4 GB/table.
- **First step.** Loop `k=304..312`, build G0 via `qmsd.triorthogonal`, call `min_dependent_columns(X_stab, 7, d_max=4)`; record d vs `needed=7`. Expect d≤3 throughout → moves the certified boundary to k=303.

### A3 · Fold p=7 m=4 into the m=2 full-span-cap crux — **RANK 3 (the only true-core closer, but it's the open problem)**

- **Idea.** The only statement that would close `[127,250]` is: *every full-rank spread S admits a full-span DGM/Leducq second-weight-class codeword of weight `≤ n/k`.* This is **exactly** the crux `D_CRUX_REDUCTION.md` isolated and left open at m=2 (verified at the p=11 extremal (21,3)-arc: true d=4 caps the line bound's 5). p=7 m=4 is strictly harder (binding word is codim-0). If the general weight-hierarchy theorem is ever proved, p=7 m=4 closes as a **free corollary**.
- **Cost / where.** **Analytic, open-ended research bet** — not a schedulable task.
- **First step.** Do **not** open it as a p=7-specific effort. Attack it once, at m=2 (smallest instance, has a verified extremal witness), as a GRM weight-hierarchy problem; treat p=7 m=4 as a downstream corollary.
- **Correction.** reformulate-P3's claim that "the m=4 distance provably reduces to 2-flats" is **false for the *punctured* distance** — min-weight *unpunctured* words are plane-supported, but the binding *punctured* words are full-span.

### A4 · Method-limitation write-up note — optional consolidation
- Convert "the bulk is hard" into a precise residual interval + a proof that the **min-weight/averaging sub-route** cannot force `d ≤ n/k` on spread S in `[127,303]` (avg config-overlap ≈ 0.009k ≈ 2.6 at k=303 ≪ the 10–15 needed). Good for a negative-result paper subsection. NB: the claimed `scratchpad/p7m4_framework.py` data does **not** exist; must be computed from scratch (hours). Scope it to the min-weight route only.

---

## Part B — CRACK IT (find a code / gather per-S evidence)

### B1 · Port the F_3 Stern/ISD finder to F_7 — **RANK 1 of the crack directions (the one genuinely new capability)**
*(merges all four lenses' top idea)*

- **Idea.** `cap_validate/stern_isd.py` is a validated, one-sided, multithreaded q-ary Stern low-weight-codeword finder — **F_3-hardcoded**. Port the numba kernel to F_7 and run it on G0 (`R=326−k` rows, `n=2401−k` cols). It works in RREF/generator form, so it **never builds `7^R`** (beats blocker #1 at any k) and finds words of **any weight** (beats #2). Every verified in-kernel codeword of weight `≤ ⌊n/k⌋` is a rigorous `d(S) ≤ n/k` witness → per-S closure, and it uniquely targets the full-span words the flat routes miss.
- **Blockers beaten.** #1 and #2 fully; #3 **only for small target weight**.
- **Cost / where.** Port ≈ **1–2 days** (mod-7 RREF, `inv[1..6]`, first-coeff normalization). Runtime set by the **true** d: minutes for `w ≤ ~9` (the `k∈[~250,303]` collapse band, **local or light lenore**); borderline core-days at `w≈11` (`k≈200`); **infeasible** at `w ≥ ~13–15` (`k ≤ ~150`). Deep sweeps are the **lenore** 32-core job.
- **First step (cheap, do before the port).** Run the *existing* overflow-free `weightcount.count_weight_d(G0, d, 7)` for `d=2,3,4` on several spread S across `[250,303]` (minutes). Since `needed(k) ≥ 7` everywhere, any weight-≤6 hit closes that S for free. **Then** port Stern to cover the `d∈[7,~9]` tail.
- **Honest framing.** Stern is **one-sided**: it can accumulate per-S closures and *hunt* a survivor, but can **never** certify `d ≥ needed`, so it cannot prove γ<1 existence and cannot produce the general no-go. "Finding nothing" near the minimum is weak evidence. Use it to (a) push the collapse-certification boundary from k=304 down toward ~250, and (b) flag any spread S it repeatedly fails to cap as an **uncertified γ<1 candidate** — never bill "no find" as closure.

---

## Dead / ruled out (one line each)

- **Recursive flat-restriction m=4 engine:** structurally omits the full-span binding words (loose: 8 vs true 3 at k=314); `structured_pe` can't even solve `RM_7(4,2)`. Loose upper bound only.
- **Extended flat/Leducq 2nd–4th class enumeration:** wrong direction — bigger supports need *more* S-concentration, which spread-S denies. Salvage = A1.
- **Bounded-unpunctured-weight enumeration for a lower bound:** correct identity, but a certified `d≥w` needs *complete* enumeration of the full-span low-weight ball (1e9–1e13 words). Certifies only the already-closed low edge.
- **GHW / Wei-duality analytic bound:** the load-bearing `max_S ↔ GHW` inequality is unproven and *is* the open problem in disguise. Circular.
- **Delsarte LP with dual-distance-294:** dual distance ≥294 is compatible with primal distance 2, so LP `max-d ∈ [21,127]`, never below `needed`. Bounds the wrong side.
- **GV / union-bound existence:** union bound is dominated by the very words that *cause* the cap; `E[killers] ≫ 1` → points to collapse, not existence.
- **Post-reduction ILP / PB-SAT lower bound:** tractable ⇒ flat-only ⇒ blind to full-span; complete ⇒ n≈2200 min-distance MIP, NP-hard and out of range.

---

## RECOMMENDED NEXT ACTION

**Ship the sure edge-tightening first, then decide.** In one local session: (A1) prove and machine-verify the coplanarity floor `d(S)≤18` and close `k∈[110,126]`; (A2) run `min_dependent_columns(d_max=4)` on `k=304..312` under `ulimit -v 4GB` and close that interval. This shrinks the certified-open window from `[110,312]` to **`[127,303]`** with zero research risk and a machine-checked certificate — a clean, publishable increment to `D_P7M4_WINDOW.md`. Only if the program wants per-S evidence across the collapse band and a genuine (slim) shot at a survivor, follow with the **B1 F_7 Stern port** (1–2 days), gated behind the free `weightcount d≤4` probe. Do **not** open A3 as a p=7 task — fold it into the existing m=2 crux.

## Honest verdict

**Marginal — worth ~1 day, not an open-ended campaign.** The window is a genuine but narrow loose end whose *hard core* `[127,~250]` is uncloseable by any current tool and is provably equivalent to the same open full-span-cap crux the Direction-D program already carries at m=2. The right terminal state is not a standalone p=7 m=4 proof but the honest status the m=2 case already holds — *edges rigorously closed, bulk capped by strong evidence (collapse + Stern), general closure deferred to the one crux that would settle all of p≤13 at once.*

---

## EXECUTION LOG (2026-07-31): A1 done (rigorous), A2 corrected

### A1 — VERIFIED (rigorous no-go, machine-checked). Closes k∈[110,126].
- Weight-21 codeword `c0` = 3 parallel lines in a coordinate 2-flat, machine-confirmed `∈ RM_7(16,4)` (orthogonal to its dual `RM_7(7,4)`).
- Lemma core verified: any 3 non-collinear points → exactly **5 of 8** plane-directions are "good" → a weight-21 codeword always meets any S in ≥3 points → **d(S) ≤ 18 for every S**.
- Explicit witness in the **real** `X_stab` at k=126: a genuine `G0^perp` codeword of weight **17 ≤ 18 ≤ n/k=18.06** ⇒ γ≥1. Crossover: `18 ≤ n/k ⟺ k ≤ 126` (fails at 127).
- **Result: certified-open window shrinks [110,312] → [127,312].**

### A2 — the dossier estimate was WRONG on both counts; NOT a cheap [304,312] closure.
Claimed: "d≤3, minutes, `min_dependent_columns(d_max=4)`." Reality on execution:
- **Distances are ~4, not ≤3.** Single random S per k: **k=312→d=4, k=311→d=4, k=310→d=4** (all γ≈1.37, γ≥1; needed d≥7 for γ<1).
- **The unbalanced MITM is memory-infeasible here:** the weight-2 table over F₇ at n≈2097 is **79 M entries (~14 GB)** and OOMs even at `d_max=5`. Only `min_dependent_columns_balanced` fits.
- **Balanced cost scales hard with dim(G0)=326−k:** ~46 s at dim 14 (k=312), but **k=304 (dim 22) exceeds 10 min** locally — the low-k end is not foreground-feasible.
- **These are random samples, not max-over-S.** Closing the band (γ≥1 for *all* S, i.e. ruling out a rare d≥7 survivor) is a **search + balanced-cert job = lenore-scale**, exactly like the mid-band. A2 is *not* a local win.

### Corrected map
- **A1 alone** shrinks the certified-open window to **[127,312]** (rigorous).
- **A2 does NOT cheaply close [304,312]** — the dossier's `[127,303]` was over-optimistic. The `[304,312]` slice needs the balanced engine + a max-over-S search (lenore), same character as the bulk. Weak local evidence (random d~4) suggests the high-k end is γ≥1, but it is unproven.
- **Certified-open window after this session: [127,312]**, with the hard core still `[127,~250]` (the m=2 full-span crux).
