# Distributed Exact Weight-Distribution Engine — Frozen Architecture

**Status:** CONTRACT for the build phases. This document is the single source of truth
for the 2-machine distributed program that certifies the **exact** minimum distance `d`
and multiplicity `A_d` (no conjecture) of the max-cap `m=7` qutrit code
`G0 = shortened RM_3(4,7)` at a size-248 cap — a `[1939, 26]` code over `F_3`
(`3^26 ≈ 2.542e12` codewords).

It merges the supplied **correctness** angle (already implemented + validated under
`dist_weightdist/`) with the **hot-core / distribution / harness** layer that the build
phases must add. Where the angle left a choice open (hot-core language, throughput
multiplier, windowed-vs-full transform on the final merge), this doc resolves it and
freezes it.

---

## 0. The computation (recap, fixed)

Given a generator `G0` of an `[n, K]` code over `F_3`:

1. **Primal weight enumerator** `A = [A_0 .. A_n]`: enumerate all `3^K` codewords
   `m -> m·G0` and histogram Hamming weights. *This is the expensive part* (`3^26`).
2. **Dual distribution** `B = [B_0 .. B_n]` of `G0^perp` via q-ary MacWilliams:
   `B_w = (1/3^K) · Σ_x A_x · K_w(x; n, 3)`, `K_w` = q-ary Krawtchouk.
   Then `d = min{ w>0 : B_w > 0 }` and `A_d = B_d`.

**Arithmetic split (the load-bearing decision):**
- Every histogram count is bounded by the message count: `A_w ≤ 3^K`. At the feasibility
  ceiling `K ≤ 29`, `3^29 = 68,630,377,364,883 ≈ 6.86e13 ≪ INT64_MAX = 9.22e18`
  (margin ≈ 134,391×; int64 stays safe to `K=39`). So **the entire enumeration +
  partial + merged histogram is int64** — no overflow, exact integer sums.
- The MacWilliams transform produces **huge** integers (`B_w` up to `~3^1913`). It is a
  **one-time** cost, negligible against `3^26`, and runs in CPython native **bignum**
  (GMP-backed `int`). **No float is ever touched** anywhere in the certified path.

---

## 1. Chosen language / tooling (RESOLVED)

| Layer | Language / tool | Rationale |
|---|---|---|
| **Hot enumeration core** | **Rust + rayon** (single static binary `dwd_core`) | Self-contained binary `scp`-able to the remote node with **no Python/toolchain dependency on the compute nodes**; rayon for intra-node thread pool; trivial trit-packing with `u64` popcount intrinsics; memory-safe so a buffer bug can't silently corrupt a partial histogram. C+OpenMP is the documented fallback (same algorithm, same `.g0`/partial formats) if a node lacks a Rust toolchain — but Rust is primary. |
| **Correctness / transform layer** | **Pure Python + native bignum** (GMP via CPython `int`) | One-time MacWilliams/Krawtchouk; reuses the trusted `qmsd.weightdist` as validation oracle; no float → a bignum can never be silently truncated. Already implemented under `dist_weightdist/` (`correctness.py`). |
| **Harness / orchestration / merge** | **Pure Python** | Work partition, checkpoint/resume, partial transfer + merge, then invokes the bignum transform. Glue, not hot. |
| **Throughput probe only** | numba-njit (`bench_kernel.py`) | A *conservative scalar floor*, not in the certified path. |

**Why not pure Python/numpy for the core:** batch matvec is `~3^26 · 1939 ≈ 5e15` ops;
measured brute int64 enumeration of `3^16 = 43M` words took `~160 s/core` in numpy →
hundreds of hours at `3^26`. The compiled core is **mandatory**.

---

## 2. Module layout under `dist_weightdist/`

```
dist_weightdist/
  __init__.py
  # --- ALREADY IMPLEMENTED + VALIDATED (correctness angle) ---
  correctness.py      # assert_int64_safe; krawtchouk_direct/_column; macwilliams (windowed+full);
                      # assert_dual_invariants; extract_d_and_Ad; merge_partials; block_checksum_ok
  io_g0.py            # read_g0 / write_g0  (.g0 v1 language-agnostic serialization)
  validate.py         # `python -m dist_weightdist.validate`  — qmsd oracle protocol (exit 0)
  bench_kernel.py     # single-core throughput floor probe
  bench/              # scratch bench data

  # --- TO BUILD (this contract) ---
  core/               # Rust crate (the hot enumeration core)
    Cargo.toml
    src/main.rs        # dwd_core binary: read .g0 (mmap), gray-code + trit-packed enumerate
    src/gray.rs        # mixed-radix base-3 gray-code message generator
    src/pack.rs        # trit-packing (2 bits/trit, 32 trits/u64) + popcount weight delta
    src/blockio.rs     # read job spec, write .partial (int64 histogram + checksum + manifest)
  core_c/              # OPTIONAL C+OpenMP fallback (same I/O contract) — build only if needed
  partition.py        # enumerate blocks; assign blocks -> nodes; emit per-node job manifests
  harness.py          # orchestrate local + (user-run) remote; spawn dwd_core; collect partials
  checkpoint.py       # checkpoint/resume state (done-blocks ledger, atomic writes)
  merge.py            # load .partial files -> merge_partials(expected_total=3^K) -> A;
                      # then correctness.extract_d_and_Ad -> {d, A_d, B}; assert_dual_invariants
  mock_cluster.py     # local two-mock-node simulation (test harness end to end, no SSH)
  formats.md          # frozen byte-level spec of .g0 / .partial / job-manifest / checkpoint
```

The `.g0` and `.partial` formats are **language-agnostic** so the Rust core and the
Python harness interoperate without serialization libraries.

---

## 3. Serialization formats (frozen)

### 3.1 `.g0` (input, already implemented in `io_g0.py`)
```
line 1:  qmsd-g0 v1 <q> <K> <n>
lines 2..K+1:  each exactly n digit-chars in {0..q-1}  (one generator row, no spaces)
```
`read_g0`/`write_g0` reduce mod q, validate shape/digits, round-trip exact. The Rust core
mmaps this and indexes each row directly. **Ingest guard:** reject a `G0` whose
`|C| = sum(A)` is not a power of 3 (rank-deficient/duplicated generator → malformed input).

### 3.2 `.partial` (per-block / per-node histogram output of `dwd_core`)
Binary, little-endian:
```
magic   : 8 bytes  "DWDP0001"
q       : u32       (=3)
K       : u32
n       : u32
nblocks : u32       number of blocks folded into this partial
checksum: u64       Σ_w hist[w]  == total messages covered  (== Σ block sizes)
hist    : (n+1) × i64   the int64 histogram A_partial[0..n]
manifest: u32 block_id list (nblocks entries) — exactly which blocks are included
```
The harness verifies `checksum == sum(hist) == Σ size(block_id)` before accepting.

### 3.3 job-manifest (harness → node) and checkpoint: see `formats.md` (built in phase 1).

---

## 4. Hot-core algorithm (`dwd_core`, Rust + rayon)

### 4.1 Gray-code message enumeration
Enumerate `F_3^K` in a **mixed-radix base-3 reflected Gray code**: consecutive messages
differ in **exactly one trit**, changing by `±1 (mod 3)`. Therefore the codeword changes
by **`± exactly one generator row`** (row `r` of `G0`), and the Hamming weight is updated
by **rescanning only the support of row `r`** — never recomputing `m·G0` from scratch.

Per step:
```
codeword[support(r)] += delta * row_r[support(r)]   (mod 3, delta ∈ {+1,+2≡-1})
update running weight by the net change in nonzeros over support(r) only
histogram[weight] += 1
```
Cost per step ≈ `|support(r)|`, not `n`.

### 4.2 Trit-packing (the speedup over the scalar floor)
Pack the length-`n` GF(3) codeword as **2 bits/trit, 32 trits per `u64`**
(`⌈n/32⌉ = ⌈1939/32⌉ = 61` words). The per-step support scan of the changed row is done
**word-parallel**: for the ≤32 trits in each touched word, add the row's trit-vector
mod 3 via the standard branch-free GF(3) word add, and derive the nonzero-count delta with
`popcount`-style masks (`is_nonzero = (lo | hi)` over the 2-bit lanes). A row of support
`s` touches `≈ s/32` words instead of `s` scalars (`~1163/32 ≈ 36` word-ops vs `1163`
scalar ops at `n=1939`, mean support_frac≈0.6).

> **Build note:** the *exact* GF(3) packed-add + nonzero-popcount lane trick is an
> implementation detail of `pack.rs`; the **contract** is only that the packed weight after
> each step equals the true Hamming weight. `validate.py` already checks the *scalar*
> reference enumeration against qmsd; the build phase MUST additionally assert
> `packed_weight == scalar_weight` on a small code before trusting the packed path.

### 4.3 Output
Each thread accumulates a **private int64 histogram** over the blocks it owns; rayon
reduces thread-private histograms into one per-node `.partial` (exact int64 vector sum).

---

## 5. int64-histogram ↔ bignum-MacWilliams split

- **int64 side (hot, distributed):** all enumeration, all partials, all merges. Enforced
  by `correctness.assert_int64_safe(K, q=3)` as a runtime contract; the Rust core uses
  `i64` and the same bound holds (proof in §0).
- **bignum side (one-time, single machine after merge):** `correctness.macwilliams`.
  Uses the q-ary Krawtchouk **three-term degree recurrence**
  `(k+1)·K_{k+1}(x) = [(q-1)(n-k)+k-q·x]·K_k(x) - (q-1)(n-k+1)·K_{k-1}(x)`,
  `K_0=1`, `K_1=(q-1)n - q·x`, with exact-integer division checked each step
  (tripwire: `num % (k+1) == 0`). Validated identical to qmsd's defining-sum formula.

**RESOLVED — which transform on the FINAL merged A:** run the **FULL** transform
(length `n+1`) on the certified merged `A`. The full path additionally asserts the global
invariant `Σ_w B_w = q^(n-K)`, giving a self-check the windowed path cannot. The
**windowed** path (`B_0..B_{d+2}`) is retained only as a fast pre-flight / for the
`extract_d_and_Ad` search guard — it **raises** (never returns a spuriously small `d`) if
true `d` exceeds the window. Final certification = full transform.

`extract_d_and_Ad(A, q, search_kmax=...)` returns
`{n, K, q, d, A_d, B, B_window, A_d_equals_Bd}`; `K` is recovered from `|C| = 3^K`, the
transform is scale-invariant under any uniform rank-deficiency multiplicity.

**Asserted dual invariants** (`assert_dual_invariants`, re-auditable standalone): every
`B_w` a non-negative integer, `B_0 = 1`, `|C|` divides `q^n`, `Σ_w B_w = q^(n-K)`.

---

## 6. Distribution / work partition

### 6.1 Blocks
Partition the `3^K` message space into **independent blocks by fixing the top `t` message
trits**: `3^t` blocks, each enumerating the remaining `K-t` trits (`3^(K-t)` messages).
For `K=26` choose `t=10` → `3^10 = 59,049` blocks of `3^16 ≈ 43M` messages each
(block granularity tunable; `t` chosen so a block runs in ~seconds-to-minutes for
fine-grained checkpointing and load balance across heterogeneous nodes).

Each block is **fully independent**: the Rust core enumerates the block's `3^(K-t)`
messages with the same gray-code/trit-packed loop, seeded at the block's base message.

### 6.2 Assignment across the two machines
`partition.py` emits a per-node job manifest (a list of block_ids). Default split is
proportional to measured per-core rate × core count (local + remote), re-balanced from the
checkpoint ledger if one node lags. The merged global enumerator is the **exact int64
vector sum** of all node partials.

### 6.3 Checksums (integrity, mandatory)
- **Per-block checksum:** a block's partial histogram sums to its message count `3^(K-t)`
  (`block_checksum_ok`).
- **Per-node checksum:** `.partial.checksum == sum(hist) == Σ block sizes`.
- **Global checksum:** `merge.py` calls
  `merge_partials(partials, n, expected_total=3^K)` which **fails loudly** on any
  lost/double-counted message. **Do NOT run the MacWilliams transform on any `A` whose
  `sum != 3^K`.**

---

## 7. Checkpoint / resume

An hours-long run must survive interruption.

- **Ledger** (`checkpoint.py`): an append-only, atomically-written record of
  `(block_id -> completed, partial_file, checksum)`. A block is "done" only after its
  `.partial` is fully written, fsync'd, and its checksum verified.
- **Resume:** on restart the harness reads the ledger, subtracts completed blocks, and
  re-dispatches only the remainder. Because blocks are independent and each writes its own
  checksum'd `.partial`, resume is exact with no double counting.
- **Atomicity:** write `.partial.tmp` then `rename()` (atomic on POSIX); the ledger entry
  is appended only after the rename succeeds. A crash mid-block leaves no ledger entry →
  the block is simply recomputed.
- **Granularity:** block size (§6.1 `t`) sets the maximum recomputation on crash (one
  block ≈ seconds–minutes).

---

## 8. Local two-mock-node test plan (`mock_cluster.py`)

**HARD CONSTRAINT honored:** never SSH / connect to `lenore_remote` or any remote host.
The two-machine topology is simulated as **two local processes** sharing the filesystem.

Test plan (all must pass before the user runs the real cross-machine job):
1. **Split correctness:** `partition.py` splits a small code's blocks across mock-node-A
   and mock-node-B (disjoint, exhaustive). Assert the union of assigned blocks == all
   `3^t` blocks, no overlap.
2. **Independent enumeration:** each mock node runs `dwd_core` on its blocks, producing
   `.partial` files. Per-block + per-node checksums verified.
3. **Transfer + merge:** harness collects both nodes' partials (simulated "transfer" = file
   copy), `merge_partials(expected_total=3^K)` succeeds, merged `A` equals the in-process
   brute enumeration of the same small code (vector equality).
4. **Checkpoint/resume:** kill mock-node-B mid-run; restart; assert the resumed run
   reproduces the identical merged `A` (and the ledger shows no recomputed-and-double-counted
   block).
5. **End-to-end:** merged `A` → full MacWilliams → `(d, A_d)` matches qmsd on the oracle
   codes (§9). Exit 0.

Recommended mock sizes: the `[[206,37,4]]_3` (K=14) and `[[667,62,4]]_3` (K=16) oracle
codes — small enough to brute-force in-process for the equality anchor, large enough to
exercise multi-block partition + merge.

---

## 9. qmsd validation protocol (`python -m dist_weightdist.validate`, exit 0)

The correctness anchor. **ALREADY PASSING** for the Python layer; the build phase extends
(c) to also run the **Rust `dwd_core`** path (not just the in-process brute enumerator) and
assert identical merged `A`.

(a) Krawtchouk recurrence (`krawtchouk_column`) == `qmsd.weightdist.krawtchouk` defining
    sum on a `(k, x, n)` grid.
(b) Full MacWilliams (`correctness.macwilliams`) == `qmsd.weightdist.macwilliams` on random
    small codes.
(c) **Oracle codes**, built via
    `qmsd.triorthogonal.build_triorthogonal_code(...)["X_stab"]`:
    - `[[206,37,4]]_3` (K=14) and `[[667,62,4]]_3` (K=16);
    - merged `A` (from the **same int64 chunked / merged-partials path** the distributed
      engine uses) == `qmsd.weightdist.weight_enumerator` EXACTLY (vector equality);
    - dual `B` == `qmsd.weightdist.macwilliams` exactly;
    - `extract_d_and_Ad` yields `d=4`, `A_d = 880` and `A_d = 3972` — matching both qmsd
      (`exact_distance_and_Ad`) and the published oracle values;
    - all dual invariants pass (`Σ_w B_w = 3^192` and `3^651`).

Because the distributed engine feeds the **identical** `extract_d_and_Ad` path, an exact
match on merged `A` guarantees an exact match on `(d, A_d)`.

---

## 10. Risks / caveats (carried from the correctness angle, binding on the build)

1. **Logical vs ALL dual codewords.** MacWilliams `B_d` counts ALL weight-`d` codewords of
   `G0^perp`; the paper's `A_d` counts only weight-`d` **logicals**. They coincide iff the
   stabilizer code has min weight `> d` (no weight-`d` stabilizers). qmsd validated this for
   all Table-3 oracles (both our cases match the published `A_d`), but it is **not provable
   from `G0` alone**. The engine reports `A_d_equals_Bd`; the real `[1939,26]` target MUST
   be checked the same way. *This is the one assumption not certified from `G0` alone.*
2. **Rank deficiency.** If `G0` is rank-deficient, `|C| = sum(A)` is a uniform multiple of
   the true count; `extract_d_and_Ad` infers `K` from `|C|=3^K` and the transform is
   scale-invariant, but a `|C|` that is **not** a power of 3 signals a malformed/duplicated
   generator — **guard at ingest** (§3.1).
3. **Windowed-search trap.** If true `d > search_kmax` the windowed path **raises** (by
   design) rather than returning a wrong small `d`. Final `d` is always confirmed by the
   **full** transform on the merged `A` (§5).
4. **Distribution integrity.** A dropped/double-counted block corrupts `A` silently except
   for the checksum. The merge MUST run with `expected_total=3^K` and every block checksum
   verified before the transform. **Never transform an `A` with `sum != 3^K`.**
5. **Packed vs scalar weight.** The trit-packed core must be asserted bit-exact against the
   scalar reference on a small code before the long run (§4.2).
6. **Krawtchouk divisibility tripwire.** Each recurrence step asserts `num % (k+1) == 0`;
   failure would indicate integer corruption (impossible in exact arithmetic) — kept as a
   tripwire.

---

## 11. Final throughput estimate — `3^26` on 52 threads (RESOLVED)

Measured single-core **floor** (numba scalar gray-step, `bench_kernel.py`, `n=1939`,
support_frac≈0.6): **≈0.89M codewords/s/core**.

| Core | Rate/core | `3^26 / (rate × 52)` | Wall time |
|---|---|---|---|
| Scalar floor (worst case) | 0.89M/s | 55,000 s | **≈15.3 h** |
| Trit-packed Rust, ~3× | 2.7M/s | ~18,100 s | **≈5.1 h** |
| Trit-packed Rust, ~6× | 5.3M/s | ~9,200 s | **≈2.6 h** |

**FROZEN ESTIMATE:** with the trit-packed Rust+rayon core the run completes in
**≈2.6–5.1 h** end-to-end; the int64 floor (if the packed speedup underdelivers) is
**≈15.3 h** — still "in hours." The bignum MacWilliams adds **seconds** (windowed
pre-flight) to **~minutes** (full final transform), negligible. **Plan budget: 6 h
target, 16 h worst-case ceiling.**

**Caveat (binding):** this assumes per-step cost ≈ row support with ≈uniform row weights.
The real punctured `RM_3(4,7)` `G0` (`n=1939, K=26`) row-weight profile sets the constant
— **re-run `bench_kernel.measure` with the actual `G0`'s mean support_frac once it is
supplied** (from the cap-construction workflow) and re-freeze the number before launch.

---

## 12. Generality

The tool reads `G0` from a `.g0` file (any `[n, K ≤ ~30]` `F_3` code). The max-cap `G0` is
supplied later. Feasibility ceiling on 52 threads ≈ `K ≤ 29` (`3^29`), where int64 still
holds with a 134,391× margin. The same pipeline (core → partition → merge → MacWilliams)
serves every code; only `bench_kernel` re-projection per-`G0` is required before a launch.
