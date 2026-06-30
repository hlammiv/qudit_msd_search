# The p=7 m=4 mid-k RM window — status and why it is hard

2026-06-30 (Direction-D item 3). The one loose end in the flat-cap no-go: is there a γ<1 code among
punctured RM₇(7,4) for k in the mid range, where the distance certifier overflows? Investigated; **partially
closed, bulk genuinely open.** This note records the structure and the validated dead-ends so they are not
re-explored.

## Structure
Code = punctured RM₇(rtilde=16, 4) on N=2401 points; dim RM(7,4)=326, so a k-puncture code is
[[2401−k, k, d]]₇ with **d ≤ d_RM = 21**. γ<1 ⟺ d > (2401−k)/k. As k grows the distance falls (more
puncturing) but the γ<1 threshold also falls, so they could cross — this window was never probed except at
high k.

## What is rigorously closed
- **Low edge, k ∈ [110,114]:** here needed-d = 21, but every nonempty puncture set meets some minimum-weight
  support (the 21-point "3 coplanar lines"), so a real codeword punctures to ≤20 < 21 ⇒ **d ≤ 20 < needed ⇒
  γ ≥ 1.** (Rigorous; min-weight codewords give a valid distance upper bound.)
- **High edge, k ≳ 304 (int64-safe, dim G0 ≤ 22):** `min_dependent_columns(d_max=3)` gives **d = 2–3** (e.g.
  k=314→3, k=320→2), far below needed-d ≈ 7 ⇒ **γ ≫ 1.** Heavy puncturing makes a large-support codeword
  almost fully punctured.

## What is open and why it is hard (the bulk, k ∈ [115, 303])
- **The min-weight shortcut fails (validated).** `d ≤ 21 − max_C|S∩supp(C)|` over the 21-point min-weight
  supports is a correct but **very loose** bound: at k=314/320 it gives d_ub=8 while the true d=2/3. The
  binding codewords are NOT the min-weight ones — they are *large* codewords whose support is almost entirely
  swallowed by S. So one cannot cheaply bound d via min-weight structure (the `maxcov` machinery in
  `scratchpad/p7m4_framework.py` is validated correct but the bound is too weak to certify the window).
- **The genuine distance is out of tool reach here.** The true d in the bulk is in the ~4–20 range. The int64
  MITM overflows for k ≤ 303 (7^(326−k) > int64). The overflow-safe (hash-encoder) MITM is `HARD_CAP`'d at
  weight 6 and is anyway infeasible at n≈2200 (a weight-6 search is ~5e11 ops). Computing d in the 7–20 range
  at n≈2200 is not feasible with any current engine.

## Assessment and recommendation
The high-k collapse to d=2–3 and the large-codeword mechanism **strongly suggest** the bulk is also capped
(γ ≥ 1), consistent with the Direction-D thesis — but it is **not certified**, and a determined search could
in principle hunt for a high-distance spread set (the verification would hit the same wall). Closing the bulk
needs either (a) a genuinely new distance algorithm reaching d≈4–20 at n≈2400 (a major, possibly infeasible
effort), or (b) **the general no-go theorem (Direction-D item 1)** — a structure-free bound that any
triorthogonal F_p code at the γ<1 density has low distance would close p=7 mid-k as a corollary, far more
cheaply and definitively than a per-k distance computation.

**Recommendation: do not pursue item 3 by direct computation** — pivot to item 1 (the general theorem),
which subsumes it. RAM note: the `min_dependent_columns` int64 MITM at dim G0 ≈ 21 transiently needs ~0.4 GB
per left table; run under `ulimit -v` (the k=305 build hit a clean MemoryError, not a crash, under the cap).
