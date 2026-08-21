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
- **`tests/`** — 161 tests anchored to the paper's published results.
- **`IMPLEMENTATION_BLUEPRINT.md`**, **`STATUS.md`** — design and honest build status.

## Highlights

- Reproduces the paper's **Table 1** (asymptotic yields), **Table 2** (smallest
  sublogarithmic codes, exact down to 18-digit block sizes), and **Table 3** search codes.
- The meet-in-the-middle routine **certifies the minimum distance of all 10 published search
  codes**, including the headline `[[519,106,5]]_5` (`d=5` in ~10 s) — distances the original
  authors' search was compute-limited to obtain. Independently adversarially verified
  (2,937 fuzz cases vs. null-space ground truth).
- **161 tests pass** (`python -m pytest -q`).

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

## Convention note

Triorthogonality and `r_max` use **`m(p-1)`** — the paper's §3.2 misprints this as `p(m-1)`;
the `m(p-1)` form is the correct, derivation-consistent one.

## Source paper

The raw arXiv paper PDF and our detailed derivation notes / undergraduate tutorial are **not
included** in this repository (the PDF is the authors' copyrighted work; the notes are kept
locally). Get the paper from [arXiv:2510.10852](https://arxiv.org/abs/2510.10852).
