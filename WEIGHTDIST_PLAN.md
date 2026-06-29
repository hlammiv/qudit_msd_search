# Implementation plan — `qmsd/weightdist.py`

An EXACT, certified MacWilliams engine. For a code whose small generator `G0` is enumerable
(`q**dim(G0)` small — the high-puncture / small-dual regime), it returns the full weight
distribution of the dual `G0^perp`, hence the certified minimum distance and `A_d`. This lifts
the `d <= 6` cap of `mindist.py` (MacWilliams cost is independent of `d`) and closes the
"exact `A_d` for large codes" gap in `STATUS.md` for that regime.

All arithmetic is in EXACT python ints (or int64 where provably non-overflowing) — never float.

---

## 0. EMPIRICAL RESOLUTION of the B_d-vs-logical-A_d question (DONE, with numbers)

Throwaway script (`scratchpad/empirical.py`) built `G0 = built["X_stab"]` for three oracle codes,
computed the dual weight distribution `B` by MacWilliams from the brute weight enumerator of `G0`
(enumerate all `q**dim` codewords), and compared `B_d` to the paper `A_d` and to the slow-but-correct
reference `qmsd.distance.A_d_logical_Z`. Results:

```
[[72,9,3]]_3 : dim(G0)=6  d=3  B_d=648   paper A_d=648   A_d_logical_Z=648   B_d==A_d? YES
              brute direct dual enum B_d=648 == MacWilliams B_d  (transform cross-checked)
[[200,43,3]]_3: dim(G0)=8  d=3  B_d=1700  paper A_d=1700  A_d_logical_Z=1700  B_d==A_d? YES
[[112,13,3]]_5: dim(G0)=7  d=3  B_d=512   paper A_d=512   A_d_logical_Z=512   B_d==A_d? YES
```

In every case **`B_d == A_d` exactly**, i.e. **the number of weight-`d` stabilizers is 0**
(`B_d - A_d_logical == 0`). All MacWilliams sanity asserts (each `B_k` a non-negative integer,
`B_0 == 1`, `sum(B) == q**n // |C|`) held, and the transform's `B_d` matched a fully independent
direct dual-codeword enumeration (648) for the 72-code.

### Why this is a THEOREM, not luck (so the engine can certify it, not just observe it)

The stabilizers are `Gp^perp`. By the standard shorten/puncture duality (`triorthogonal.py` docstring,
NOTES Thm 3): `G0^perp = puncture(RM(rtilde,m), S)` (the big logical-carrying code) and
`Gp^perp = shorten(RM(rtilde,m), S)`. **Shortening never decreases minimum distance**, so every
nonzero stabilizer has weight

```
    wt >= d_min(Gp^perp) >= d_RM(rtilde, m)       (closed form: qmsd.reedmuller.d_rm)
```

Therefore **if the logical distance `d < d_RM(rtilde, m)`, there are zero weight-`d` stabilizers, and
`A_d_logical == B_d` EXACTLY and CERTIFIABLY** (no enumeration needed). Checked for all six oracle
codes — `d` is always strictly below `d_RM(rtilde,m)`:

```
[[20,5,2]]_5 : rtilde=5 d_RM(5,2,5)=4 > 2     [[72,9,3]]_3 : rtilde=5 d_RM(5,4,3)=6 > 3
[[112,13,3]]_5: rtilde=8 d_RM(8,3,5)=5 > 3    [[200,43,3]]_3: rtilde=6 d_RM(6,5,3)=9 > 3
[[206,37,4]]_3: rtilde=6 d_RM(6,5,3)=9 > 4    [[667,62,4]]_3: rtilde=8 d_RM(8,6,3)=9 > 4
```

### Decision for `exact_distance_and_Ad`

* `B_d = B[distance]` is always computed (cheap, certified, multiplicity-invariant).
* `A_d_logical` defaults to `B_d`. This is EXACT whenever there are no weight-`d` stabilizers, which is
  guaranteed by the `d < d_RM(rtilde,m)` test — true for the entire paper regime.
* `A_d_equals_Bd` records whether that equality is certified. The bare engine (which receives only `G0`,
  not `Gp`/`rtilde`) sets `A_d_equals_Bd = True` under the documented no-weight-`d`-stabilizer regime; the
  **`codes.py` integration** (which has the full `built` dict + `m, r` → `rtilde`) performs the actual
  `d < d_RM(rtilde,m)` certification and, when feasible (`C(n,d)` within budget), cross-checks against
  `distance.A_d_logical_Z`. The engine must NEVER silently report `B_d` as `A_d` when they could differ:
  if a caller cannot establish `d < d_RM(rtilde,m)`, it must treat `A_d_logical = B_d` as an UPPER BOUND
  and fall back to `A_d_logical_Z` (the regime where this matters — `d == d_RM` — does not occur in the
  paper codes but is handled honestly).

---

## 1. Module layout (`qmsd/weightdist.py`)

```
module docstring          # the math (Krawtchouk, MacWilliams), the regime, the B_d==A_d theorem above
imports                   # numpy, math.comb, functools.lru_cache, galois only via reedmuller helpers

_krawtchouk_cached(...)   # @lru_cache inner, exact-int
def krawtchouk(k,x,n,q)   # public, thin wrapper over the cache

def macwilliams(A, q)     # exact-int transform + the three asserts

def _full_rank_rows(G,q)  # row-reduce G over F_q, drop zero rows -> (rank x n) int64  (engine-internal)
def weight_enumerator(G, q, max_words=5_000_000)   # chunked numpy histogram [A_0..A_n]
def dual_weight_distribution(G0, q, max_words=...)  # macwilliams(weight_enumerator(...))
def exact_distance_and_Ad(G0, q, max_words=...)     # the engine dict
```

No new heavy package-level re-exports (keep `qmsd/__init__.py` import-light, per its docstring).

---

## 2. Build order

1. **`krawtchouk(k, x, n, q)`** — `sum_{j=0}^{k} (-1)**j (q-1)**(k-j) C(x,j) C(n-x,k-j)` with
   `math.comb`, pure exact int. Memoize on `(k,x,n,q)` (`lru_cache`); the transform calls it
   `(n+1)**2` times so caching matters. Smoke vs hand values: `K_0 = 1`; `K_1(x;n,q) = (q-1)n - q x`.
2. **`macwilliams(A, q)`** — `n = len(A)-1`, `C = sum(A)` (python int). For each `k`,
   `s = sum(A[x]*krawtchouk(k,x,n,q) for x in range(n+1) if A[x])`; `assert s % C == 0`;
   `B_k = s // C`. Then `assert B[0]==1`, `assert all(b>=0)`, `assert (q**n) % C == 0 and
   sum(B) == q**n // C`. (These hold provided `C == q**dim` — see the full-rank note in §3.)
3. **`weight_enumerator(G, q, max_words)`** — chunked numpy enumeration (next section).
4. **`dual_weight_distribution`** — one-liner: `macwilliams(weight_enumerator(G0,q,max_words), q)`.
5. **`exact_distance_and_Ad`** — feasibility guard, then distance / `B_d` / `A_d_logical` per §0.
6. **`codes.py` integration** — opportunistic exact-`A_d`/distance path when `dim(G0)` is small
   (§5), gated so it never changes existing behavior or breaks tests.

---

## 3. `weight_enumerator` — exact, chunked, numpy, never materialize all codewords

Contract: `weight_enumerator(G, q, max_words=5_000_000) -> [A_0..A_n]`, the exact weight enumerator of
the code **spanned by `G` (`r x n`)**, by enumerating all `q**r` messages in chunks and accumulating a
length-`(n+1)` histogram. `raise ValueError` if `q**r > max_words`.

Algorithm:
* `Gi = (np.asarray(G).astype(int64)) % q`; `r, n = Gi.shape`. `total = q**r`; if `total > max_words`
  raise `ValueError` (with the `q**r` / `max_words` numbers in the message).
* Iterate message indices `0..total-1` in chunks of `CHUNK = min(total, ~200_000)`:
  - For a chunk of message ints `idx` (an int64 array), extract the `r` base-`q` digit columns by
    repeated `idx // q**j % q` → coefficient matrix `M` of shape `(chunk, r)`, entries in `[0,q)`.
  - `words = (M @ Gi) % q` → `(chunk, n)`. **Overflow-safe in int64**: each entry before `% q` is
    `<= r*(q-1)**2` (e.g. `274*16 ≈ 4400`), far under int64. (For exotic huge `r` the plan allows an
    object-dtype fallback, but the `max_words` guard makes this irrelevant in practice.)
  - `w = np.count_nonzero(words, axis=1)` → per-row Hamming weights.
  - `A += np.bincount(w, minlength=n+1)` accumulated into a python-int list at the end.
* Return `A` as a list of python ints. Memory is `O(CHUNK * n)` only.

**Multiplicity / full-rank note (a real subtlety — documented, not glossed):** enumerating `q**r`
messages counts each *distinct* codeword `q**(r-rank)` times. For the MacWilliams `B` this multiplicity
`λ = q**(r-rank)` cancels (`B_k = (1/λ|C|) sum λ A_x K_k = (1/|C|) sum A_x K_k`), so `B` is correct
regardless. BUT the `sum(B) == q**n // sum(A)` assert needs `sum(A) == q**dim`, i.e. `r == rank`.
So **the engine (`dual_weight_distribution`/`exact_distance_and_Ad`) row-reduces `G0` to full row rank
first** (`_full_rank_rows`, via `galois` `row_reduce` + drop zero rows, exactly as `distance._code_basis`
does) before calling `weight_enumerator`. `weight_enumerator` itself stays pure (no row-reduce) but its
docstring states it assumes/expects independent rows for the downstream sanity asserts to be exact.
(The oracle `X_stab` matrices have full row rank = the quoted `dim(G0)`, confirmed empirically.)

---

## 4. `exact_distance_and_Ad(G0, q, max_words)` — the engine

Returns `{feasible, distance, weight_dist, B_d, A_d_logical, A_d_equals_Bd}`.

* `R = _full_rank_rows(G0, q)`; `dim = R.shape[0]`.
* **Feasibility:** if `q**dim > max_words` → return
  `{feasible: False, distance: None, weight_dist: None, B_d: None, A_d_logical: None,
    A_d_equals_Bd: None}` (no enumeration attempted).
* Else `B = dual_weight_distribution(R, q, max_words)` (asserts fire inside `macwilliams`).
* `distance = min(w for w in range(1, n+1) if B[w] > 0)`. (`B[0] == 1` asserted; a nonzero dual exists
  because `dim < n` in the regime — if `B[1:]` were all zero the dual is trivial; guard: if no `w>0`
  has `B[w]>0`, return `distance=None` honestly.)
* `B_d = B[distance]`.
* `A_d_logical = B_d`; `A_d_equals_Bd = True` — the documented no-weight-`d`-stabilizer regime
  (§0). The function docstring states the certification condition `d < d_RM(rtilde,m)` and that the
  `codes.py` layer enforces it; `weight_dist = B`.

---

## 5. `codes.py` integration (additive, non-breaking)

Add an opportunistic exact path that does NOT change any existing default behavior (so the 165 tests
stay green). Concretely, a new helper (e.g. `exact_Ad_via_macwilliams(built, p, m, r, max_words=...)`)
or an optional branch inside `code_from_puncture`:

* Compute `dim = rank(G0)`. If `p**dim <= max_words` (small-dual regime), call
  `weightdist.exact_distance_and_Ad(G0, p)`; obtain `distance`, `B_d`.
* Certify `A_d_equals_Bd` here, where `rtilde = reedmuller.r_tilde(r, m, p)` is available:
  `cert = distance < reedmuller.d_rm(rtilde, m, p)`. If `cert`, set `A_d = B_d` (EXACT, certified);
  set `Code.d` to `distance` even when `distance > 6` (this is the cap-lifting win over `mindist`).
* If `not cert` (would only happen outside the paper regime), fall back to
  `distance.A_d_logical_Z(built, p, distance)` when `C(n,distance)` is within its budget; otherwise
  report `B_d` as an explicit upper bound and leave `A_d` unset/flagged — never a silent guess.
* When `p**dim > max_words`, leave the current behavior untouched (MITM distance + `A_d_logical_Z`
  fallback already in `code_from_puncture`).

This routes the high-`k` / small-dual codes (exactly the regime DISTANCE_ALGORITHM_RESEARCH.md calls
the "winner's corner") through MacWilliams for certified `d` AND `A_d`, including `d > 6`.

---

## 6. Test plan — `tests/test_weightdist.py`

1. **Krawtchouk unit values.** `krawtchouk(0,x,n,q) == 1`; `krawtchouk(1,x,n,q) == (q-1)*n - q*x`;
   `krawtchouk(k,0,n,q) == comb(n,k)*(q-1)**k`; symmetry/recurrence spot checks. A couple of fully
   hand-computed entries for `q in {2,3,5}`.
2. **MacWilliams involution + self-dual sanity.** For random small generators `G` (`r x n`,
   `q in {2,3,5}`, `n<=10`, `r<=5`): compute `A = weight_enumerator(G,q)` (full-rank reduced),
   `B = macwilliams(A,q)`, and `macwilliams(B', q)` should return `A` scaled appropriately
   (`MacWilliams` of `MacWilliams` recovers the primal up to `|C|` bookkeeping). At minimum assert
   the three invariants and `B_0==1`.
3. **`weight_enumerator` brute cross-check.** For random small codes, compare the chunked-numpy
   histogram against an independent brute enumeration (the `distance.py`-style explicit codeword loop)
   — exact list equality. Include a rank-deficient `G` to exercise the multiplicity behavior, and a
   case with `q**r` straddling a small `max_words` to assert `ValueError` fires.
4. **Dual distribution vs `mindist` distance.** For several random small `G0`: the smallest `w>0`
   with `dual_weight_distribution(G0,q)[w] > 0` equals `mindist.min_dependent_columns(G0, q)` (the
   independently-verified MITM). This is the MITM distance cross-check, and it also exercises `d`
   beyond the MITM's own `d<=6` cap on at least one constructed code (MacWilliams has no cap).
5. **ORACLE ground truth (the headline test).** For the six small-dual oracle codes
   (`[[20,5,2]]_5, [[72,9,3]]_3, [[112,13,3]]_5, [[200,43,3]]_3, [[206,37,4]]_3, [[667,62,4]]_3`):
   build `G0 = built["X_stab"]`, run `exact_distance_and_Ad(G0, p)` and assert
   `distance == paper d` and `B_d == A_d_logical == paper A_d` (from `tests.ground_truth.TABLE3_CODES`).
   `[[667,62,4]]_3` has `dim(G0)=16` → `3**16 ≈ 4.3e7` messages; keep it (chunked, ~seconds) or mark
   `slow` if CI time matters.
6. **`B_d == A_d_logical_Z` consistency.** For at least the `dim<=8` oracle codes, assert the engine's
   `A_d_logical` equals `distance.A_d_logical_Z(built, p, distance)` (the slow reference) AND the
   certification condition `distance < reedmuller.d_rm(rtilde, m, p)` holds — pinning the §0 theorem.
7. **Feasibility guard.** A code with `q**dim(G0) > max_words` (small `max_words`) returns
   `feasible=False` with all-`None` fields and performs no enumeration.
8. **`codes.py` integration test.** `code_from_puncture` (or the new helper) reproduces `d` and `A_d`
   for the small-dual oracle codes via the MacWilliams path, matching `TABLE3_CODES`.

Run: `python -m pytest tests/test_weightdist.py -q`, then the full suite `python -m pytest -q`
(must stay green; do not modify files other than `qmsd/weightdist.py`, the additive `codes.py` hook,
and the new test file).

---

## 7. Exact aggressive-check assertions (must all be present)

In `macwilliams` (the core certifier):
* `assert s % C == 0` for every `k` (each `B_k` is an exact integer) — fail loud otherwise.
* `assert B[0] == 1`.
* `assert all(b >= 0 for b in B)`.
* `assert (q**n) % C == 0`.
* `assert sum(B) == q**n // C`  (`== q**(n-dim)` when `C == q**dim`).

In `weight_enumerator`:
* `assert q**r <= max_words` else `ValueError` (the feasibility contract).
* `assert A[0] >= 1` (zero codeword present) and `assert sum(A) == q**r` (every message counted).
* int64 non-overflow is structural (`<= r*(q-1)**2`); assert `r*(q-1)**2 < 2**62` as a guard, else
  fall back to object dtype.

In `exact_distance_and_Ad`:
* On `feasible=True`: `assert weight_dist[0] == 1`; `assert distance is None or 1 <= distance <= n`;
  `assert B_d == weight_dist[distance]`.
* `assert A_d_logical == B_d` together with `A_d_equals_Bd is True` in the certified regime; the
  `codes.py` caller additionally `assert distance < d_rm(rtilde, m, p)` before trusting `A_d = B_d`,
  and otherwise routes to `A_d_logical_Z` — guaranteeing `B_d` is never silently mislabeled as `A_d`.

---

## 8. Files

* CREATE `qmsd/weightdist.py` (the engine).
* CREATE `tests/test_weightdist.py` (the test plan above).
* EDIT (additive, non-breaking) `qmsd/codes.py` — opportunistic MacWilliams exact path for the
  small-dual regime (§5).
* Do NOT touch any other module; keep all 165 existing tests green.
