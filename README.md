# qudit_msd_search

Discovering low-overhead **quantum error-correcting codes for qudit magic state
distillation (MSD)** — a reimplementation and extension of the search in
[arXiv:2510.10852](https://arxiv.org/abs/2510.10852), *"Sublogarithmic Distillation in all
Prime Dimensions using Punctured Reed-Muller Codes"* (Tanay Saha & Shiroman Prakash, 2025).

Given a prime qudit dimension `p`, the `qmsd` package builds triorthogonal codes by
puncturing maximal-degree Reed-Muller codes over `F_p` and scores them by the yield
parameter `gamma = log(n/k)/log(d)` (overhead exponent; `gamma < 1` is *sublogarithmic*)
and by single-round distillation cost `C`.

## What's here

- **`qmsd/`** — the discovery package. Two engines: an exact integer **analytic** engine for
  the Manhattan family (scales to astronomically large codes) and an explicit **finite-field**
  engine for searching arbitrary puncture sets. Includes `mindist.py`, a meet-in-the-middle
  exact minimum-distance routine, and `data/puncture_locations.json`, the machine-validated
  puncture columns of the paper's 10 search codes (the correctness oracle).
- **Distillation performance** (following
  [arXiv:2610.08012](https://arxiv.org/abs/2610.08012), Prakash, *Cyclic triorthogonal codes in
  prime dimension*): `performance.py` (exact output error `eps_out(eps)`, acceptance
  probability, leading order `c·eps^d` and single-round threshold, per output and per block),
  `protocol.py` (multi-round protocol cost and rigorous bounds when exact enumeration is out of
  reach), `cyclic.py` (the paper's cyclic constructions A/B), `mc_distill.py` (Monte Carlo
  cross-check), and `dist_weightdist/core_c/dwq_coset.c` (multithreaded C enumerator for large
  codes). See **`PERFORMANCE.md`** for the model, validation and results.
- **Doubled codes and weighted supports** (following
  [arXiv:2610.29457](https://arxiv.org/abs/2610.29457), Prakash & Ramya, *Short qudit
  triorthogonal codes from doubling*): `doubling.py` (P-codes, Theorem 4 doubling, the paper's
  Table 3 codes and the Prakash–Saha [[6m+2,3m−2,2]]₃ family, all in the standard protocol
  library), `objective.py` (search objective: a code's protocol-cost saving over the standard
  codes, with a fast leading-order surrogate) and `weighted.py` (weighted Reed–Muller supports,
  i.e. shorter triorthogonal spaces with charge coefficients). Notes:
  `literature/2610.29457_NOTES.md`; results: **`WEIGHTED_SUPPORTS.md`** (the first γ < 1
  search code at p = 7, its certificate and the class-block lemma behind it).
- **`SMALL_CODES.md`** — small effective codes at every prime p ≤ 97: the full Reed–Solomon family
  with exact thresholds and coefficients (closed-form MDS engine, `qmsd.rs_family`), plus the small
  search and literature codes.
- **`CATALOG.md`** — the catalog record schema (v2): every code carries a rebuild recipe, a
  structured distance certificate (bounds, methods, witness), performance pointers and geometry;
  `qmsd catalog show <id>` prints them. `CERTIFY_P7.md` tracks the running certification queue.
- **`tests/`** — 400 tests anchored to the published results of the papers.
- **`IMPLEMENTATION_BLUEPRINT.md`**, **`STATUS.md`** — design and honest build status.

## Highlights

- Reproduces the paper's **Table 1** (asymptotic yields), **Table 2** (smallest
  sublogarithmic codes, exact down to 18-digit block sizes), and **Table 3** search codes.
- The meet-in-the-middle routine **certifies the minimum distance of all 10 published search
  codes**, including the headline `[[519,106,5]]_5` (`d=5` in ~10 s) — distances the original
  authors' search was compute-limited to obtain. Independently adversarially verified
  (2,937 fuzz cases vs. null-space ground truth).
- **The first γ < 1 magic-state distillation code at p = 7: [[1756,260,7]]₇, γ = 0.982**, a
  weighted Reed–Muller support (sphere charge on F₇⁴) with outputs on one charge class, d = 7
  certified by a three-part proof plus an explicit weight-7 word (`WEIGHTED_SUPPORTS.md`). The
  best published p = 7 codes have γ ≥ 1.6; the sibling [[1746,270,6]]₇ cuts distillation cost
  ~4× below the cyclic and doubled p = 7 codes in the leading-order surrogate.
- Numba meet-in-the-middle kernels certify distances to d = 7 at n ~ 2000 (bucket-indexed, division-free
  right scans: 6-7x the probe rate of a binary-search join at identical results, `scratchpad/bench_mitm/`).
- Reproduces arXiv:2610.08012 **Table 2** (leading order, threshold, acceptance) to the printed
  digits, and gives exact thresholds or rigorous threshold lower bounds for the catalog codes
  (`PERFORMANCE.md`).
- **tests: see `python -m pytest -q`** (400+ collected, 7 skipped) (`python -m pytest -q`).

## Quickstart

Requires Python ≥ 3.10 with `numpy`, `galois`, `scipy`. From the repo root:

```bash
python -m qmsd search --p 5 --m 4 --trials 5000      # search for low-gamma codes
python -m qmsd reconstruct --label "[[519,106,5]]_5" # rebuild + certify a published code
python -m qmsd asymptotic --p 5                      # asymptotic optimal yield gamma_0(p)
python -m pytest -q                                  # run the test suite
```

## Add new search results to the explorer

Export the displayed search candidates as structured JSON, then import them into the
local explorer catalog:

```bash
python -m qmsd search --p 3 --m 5 --sampler flat_spread --target-k 13 \
  --trials 300 --jobs -1 --output runs/p3-m5-k13.json
python -m qmsd catalog import runs/p3-m5-k13.json
streamlit run app.py
```

Imports are validated, deduplicated under stable artifact IDs, and written to
`qmsd/data/catalog/`. Certified, full-rank search results are inferred as `confirmed`;
uncertified results are imported as `candidate` and cannot be promoted to `confirmed`
by a command-line flag. The running explorer notices catalog changes on its next rerun.

Use `--status candidate` to deliberately keep a certified result provisional, or
`--catalog-dir PATH` (equivalently `QMSD_CATALOG_DIR`) for a separate catalog.

## Distillation performance and protocol cost

Beyond `gamma`, each catalog code can be scored by how it actually distills under twirled
depolarizing noise (paper convention `eps = (p-1)/p · delta`):

```bash
python -m qmsd perf <artifact_id>                    # leading order + exact threshold for one code
python -m qmsd perf --all --dry-run                  # which codes are exact locally (tier 1)
python -m qmsd assess <artifact_id>                  # threshold, P_acc, multi-round cost vs standard codes
python -m qmsd value <artifact_id> [--surrogate]     # search objective: cost saving over the standard library
python -m qmsd assess --all --store --skip-existing  # assess every new catalog code; shown in the explorer
python -m qmsd protocol --p 5 --eps0 1e-2 --target 1e-15  # cheapest multi-round protocol for p
```

`assess` prints the cost (raw magic states per output) of the cheapest protocol with and
without the code, compared with the standard codes: Reed–Solomon, Reed–Muller, cyclic, doubled,
and (p = 3) Prakash–Saha. Against the Prakash–Saha codes our qutrit codes save 1.0–1.6×.
Small codes get exact numbers. Larger codes get rigorous bounds, with `≥` marking a bound,
from an enumerable subcode (Lemma 17 of arXiv:2610.08012). Results are stored under
`qmsd/data/performance/` and appear in the explorer's Distillation tab and the "best saving"
column.

## Convention note

Triorthogonality and `r_max` use **`m(p-1)`** — the paper's §3.2 misprints this as `p(m-1)`;
the `m(p-1)` form is the correct, derivation-consistent one.

## Source paper

The raw arXiv paper PDF and our detailed derivation notes / undergraduate tutorial are **not
included** in this repository (the PDF is the authors' copyrighted work; the notes are kept
locally). Get the paper from [arXiv:2510.10852](https://arxiv.org/abs/2510.10852).
