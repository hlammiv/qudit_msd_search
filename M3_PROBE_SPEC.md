# Probe spec: can m=3 give small γ<1 codes at p=11/13? (de-risk before the full build)

**Status:** EXECUTED 2026-06-29 → **VERDICT: NO-GO.** A *bounded, verifiable* experiment to decide whether
the full m=3 structured-enumerator build (the project's biggest piece) is worth committing to.

## RESULT (2026-06-29) — NO-GO, brute-confirmed
Built + validated `qmsd/structured_m3.py` (line+plane enumerator): validation 2a (m=3 lines restricted to a
plane == structured_pe m=2) PASSES; validation 2b (d_struct == certified true distance) PASSES for p=5/7
(MacWilliams cross-checked) AND at p=11 (structured `d_planes=3` matched brute `min_dependent_columns`=3).
**At p=11 m=3 in the γ<1 window the TRUE distance is 2–3** (k=210→d=3 γ=1.525; k=212→3; k=214→2; k=216→2),
far below the d=6 needed. Lines give d_lines=6 (max_line=5, γ≈0.94) but **the plane (2-flat) codewords cap
d at 2–3**: at the γ<1 density each plane (p²=121 pts) catches ~k/11≈19 punctures, so plane codewords
puncture to weight 2–3 regardless of flat-spread. Same mechanism as the m=2 2D-cap, more severe.
⇒ **p=11/13 are plane-capped at m=3 too; the full Kasami–Tokura build is NOT justified. Banked.** The
line-spread insight works only where line codewords are the true minimum (p≥17 m=2). `qmsd/structured_m3.py`
+ tests are a kept, validated deliverable.

## 0. Why this probe exists
- m=2 is **2D-capped** for p=11/13: even with a max_line≤3 puncture set (ILP-confirmed to exist), a 2D
  full-span codeword punctures to ≤4, so d≤4 ⇒ γ>1. The line-spread cannot cross at m=2. (Settled 2026-06-29.)
- m=3 has more distance room (block p³, larger d_RM), and the paper's p=11/13 γ<1 codes are analytic/huge —
  so *small* m=3 γ<1 codes would be new.
- **But at m=3 the brute MITM is infeasible** (the γ<1 window has n~10³; weight-3 enumeration is
  `C(n,3)·(p−1)² ≈ 1.6e10` ⇒ ~500 GB even rebalanced). So a structured (flat) enumerator is the *only* route,
  and there is **no brute oracle** to validate large cases. That is the risk this probe measures.

## 1. The question the probe answers (go/no-go)
**At p=11, m=3, in the γ<1 window, can a flat-spread puncture set make the LINE+PLANE structured distance
large enough to clear γ<1 — or do the plane (2-flat) codewords cap it (the m=3 analog of the m=2 2D cap)?**

- If even the optimistic line+plane distance **cannot clear γ<1** → **NO-GO** (higher Kasami–Tokura classes
  only lower it further; don't build the full thing).
- If it **clears γ<1** → **PROMISING**; commit to the full build (all weight classes + rigorous certification).

This is an *upper bound* on what's achievable; a clean negative here kills the project cheaply, a positive
justifies the expensive build.

## 2. The math (m=3, the test case)
For p=11, m=3: `r_max = ⌊(3·10−1)/3⌋ = 9`, `rtilde = 3·10 − 9 − 1 = 20 = 2·10 + 0` ⇒ `a=2, b=0`,
`d_RM = (p−b)·p^(m−a−1) = 11·11⁰ = 11`. The minimum-weight RM₁₁(20,3) codewords are **line-supported**
(1-flats, p points each) — same structure as `structured_pe` but in AG(3,p). The next class up is
**plane-supported** (2-flats, p² points), the m=3 analog of the m=2 "2D" codewords that did the capping.

So the probe needs the **flat distance from 1-flats AND 2-flats**:
- `d_lines  = d_RM − max_{line ℓ}  |S∩ℓ|`   (1-flat contribution — generalize `structured_pe.enumerate_lines`/
  `line_punctured_distance` from AG(2,p) to AG(3,p))
- `d_planes = (plane min-weight) − (plane-puncturing)` — enumerate the 2-flats of AG(3,p) and the RM codewords
  supported on each (a plane restricts to RM_p(rtilde−(p−1), 2), i.e. the m=2 problem we already solve), then
  the punctured weight `|supp(c)\S|` over plane-supported `c`.
- `d_struct = min(d_lines, d_planes)` is the probe's distance estimate (an **upper bound** on the true d).

**References:** Delsarte–Goethals–MacWilliams IC16(1970) (min-weight = affine flats); Kasami–Tokura IC30(1976)
(weights < 2.5·d_min); Leducq arXiv:1001.2554. Reuse the decomposition-by-minimal-affine-span idea already in
`qmsd/structured_ad.py` (p=3) and `qmsd/structured_pe.py` (p=17 m=2, lines).

## 3. Tasks (in order; stop at any NO-GO)
1. **Build** `qmsd/structured_m3.py` (additive, new file): `flat_distance(p, 3, r, puncture_columns)` returning
   `{d_lines, d_planes, d_struct, witness}` — enumerate 1-flats and 2-flats of AG(3,p), compute punctured
   weights, no `(p−1)^a·C(n,a)` materialization. Restrict-to-plane reuses the existing m=2 line logic.
2. **Validate** (decisive — this is the whole point, since no brute oracle at scale):
   - **(a)** The m=3 line logic reproduces `structured_pe` exactly when restricted to a plane (the m=2 case).
   - **(b)** On a SMALL p=11 m=3 code with `dim(G0) ≤ 8` (so `11^dim(G0) ≤ 2.4e8` MacWilliams is feasible,
     `qmsd.weightdist.exact_distance_and_Ad`), `d_struct` must equal the true distance **exactly**. If it
     doesn't, the enumerator is incomplete → fix or report the gap before going further.
3. **Flat-spread sampler:** minimize `max_flat |S∩flat|` over BOTH lines and planes (ILP via `scipy.optimize.milp`
   — the AG(3,p) incidence matrix; or annealing). Confirm what min flat-intersection is reachable at the
   γ<1 densities.
4. **The probe run:** sweep the γ<1 window (high k, `dim(G0)` small), flat-spread sample, compute `d_struct`,
   and report the best γ. Decide GO/NO-GO per §1.

## 4. Success / failure criteria
- **Validation gate:** if step 2(b) fails (d_struct ≠ MacWilliams d on the small case), STOP — the enumerator
  is wrong/incomplete; report and do not trust the probe numbers.
- **GO:** validation passes AND some flat-spread p=11 m=3 code has `d_struct` giving γ<1 (with the caveat that
  the full classification must still confirm no higher-weight class drops it). ⇒ write the full-build spec.
- **NO-GO:** validation passes but no flat-spread code clears γ<1 even with the optimistic line+plane distance.
  ⇒ p=11/13 are 2D/plane-capped at m=3 too; bank the whole effort. Clean negative.

## 5. Scope boundaries (keep it bounded)
- **IN:** 1-flats and 2-flats only; the p=11 m=3 test; one MacWilliams validation; the GO/NO-GO call.
- **DEFERRED to the full build (only if GO):** the complete Kasami–Tokura classification (line/plane *unions*
  and the intermediate weight classes up to 2.5·d_min); rigorous d-certification with no brute fallback; the
  flat-spread A_d; p=13, p=7-at-m≥3.
- **Do NOT** attempt the brute MITM at m=3 (infeasible — §0) or materialize any `C(p³,a)` array.

## 6. Deliverable + estimate
- Deliverable: `qmsd/structured_m3.py` + `tests/test_structured_m3.py` + a one-paragraph GO/NO-GO verdict with
  the validation result and the best p=11 m=3 γ.
- Estimate: ~1 day / one multi-agent workflow. Risk concentrated in step 2(b) (validation) and step 4 (whether
  the planes cap it). This is the cheap experiment that decides the expensive build.
