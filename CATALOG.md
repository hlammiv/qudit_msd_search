# The code catalog (schema `qmsd.catalog-record.v2`)

*2026-10-29. Records: `qmsd/data/catalog/*.json`, loaded by `qmsd.results.load_result_catalog()`;
built-in records (paper tables, analytic Manhattan family, top-level artifacts) are generated in
code with the same fields. CLI: `qmsd catalog import | migrate | show`.*

## Why v2

v1 records stored a code as "p, n, k, d, a puncture list and a sentence of evidence". That could
not rebuild the weighted-support codes (the paper's p = 7 headline is not a punctured RM code), and
"d = 7 exact with witness", "d >= 7 with the upper bound open" and "closed-form MDS distance" were
indistinguishable to a program. v2 keeps every v1 field and adds four structured blocks. The
promise of the catalog is: **every record can be rebuilt from its `construction`, and its claim is
exactly what its `certificate` says.** `tests/test_catalog_v2.py` checks both for every record.

## Record layout

| field | meaning |
|---|---|
| `schema` | `qmsd.catalog-record.v2` (v1 files are upgraded on load; `qmsd catalog migrate` rewrites them) |
| `artifact_id` | stable id, hash of the construction identity (`imported-p{p}-n{n}-k{k}-{hash}` for imports) |
| `label`, `p`, `m`, `r`, `w`, `n`, `k`, `d`, `A_d`, `gamma_value` | summary parameters; `d` is the certificate's upper bound when the bounds meet, else the lower bound |
| `puncture_columns` | 1-indexed puncture columns (punctured RM / RS only; kept for v1 tools) |
| `family`, `provenance`, `source`, `note` | human-readable origin |
| `status` | `confirmed` / `partial` / `candidate` / `refuted`, derived from the certificate (below) |
| `distance_evidence`, `Ad_evidence` | one-sentence summaries (v1 compatibility; the certificate is authoritative) |
| `construction` | rebuild recipe (next section) |
| `certificate` | distance bounds, methods, witness, stages |
| `performance` | which performance data exist |
| `geometry` | optional search-family descriptors |

### `construction`

```
{"family": "punctured_rm" | "reed_solomon" | "weighted_rm" | "manhattan" | "none",
 "p": 7, "m": 4, "r": 7,
 "puncture_columns_1indexed": [...],            # punctured_rm, reed_solomon (m = 1)
 "w": 20,                                       # manhattan: puncture all points of Manhattan weight <= w
 "charge": {"kind": "quadric", "diag": [1,1,1,1], "lin": [0,0,0,0], "const": 0},   # weighted_rm
 "outputs_in_X": [...], "class_unit": 1}        # weighted_rm: 0-based indices into the support X
```

`qmsd.performance.record_matrices(record)` dispatches on `family` and returns (G0, H1): the
shortened generator (X-stabilisers, the acceptance test) and the canonical logical rows. Weighted
supports go through `qmsd.weighted.weighted_support` + `weighted_code` with the stored charge
(`quadric(p, m, diag, lin, const)`) and output set; the sweep and certification scripts write
exactly these fields. `family: none` marks records that cannot be rebuilt (the refuted cap code,
the astronomically large Manhattan records are rebuilt analytically instead).

### `certificate`

```
{"lower_bound": 7, "upper_bound": 7, "exact": true,
 "lower_method": "no weight <= 5 (MITM + parallel weight-5 scan); no weight-6 block word (3+3); no mixed weight-6 word (line lemma)",
 "upper_method": "explicit weight-7 word (3+4 block scan)",
 "witness": {"columns": [523, 1045, ...], "coefficients": [3, 3, 1, ...], "weight": 7},   # G0 columns, 0-based
 "stages": [{"stage": "block6", "result": "...", "seconds": 4375}, ...],
 "machine": "lenore, 8 threads nice 19", "date": "2026-10-29", "refuted": false}
```

* `lower_bound` is what has been *proved* about the minimum distance; `upper_bound` is what an
  explicit word (or a closed form) shows. `exact` iff they meet.
* The status is derived (`qmsd.catalog.status_from_certificate`): exact -> `confirmed`; a rigorous
  lower bound with the exact value open, or a closed-form distance not independently re-certified
  -> `partial`; no rigorous bound -> `candidate`; `refuted` when a later check found a lighter word.
* A witness, when stored, is verified against the rebuilt G0 by the test suite (`G0 w = 0 mod p`
  and the stated weight). Records imported from the search engine carry the engine's MITM
  certificate (exact, no stored witness).

### `performance`

```
{"exact_curve": false, "rigorous_bounds": false,
 "census": {"A_d": 24, "L_d_per_output": 1.33, "c_eps": 2.9e-5, "source": "runs/weighted/...census.json"} | null,
 "acceptance_mc": {"source": "runs/certify/pacc_mc.json", "P_acc": {"0.001": 0.174, ...}} | null,
 "notes": "exact enumeration infeasible (dim A = 66)"}
```

Pointers, not copies: the exact curves and thresholds live in `qmsd/data/performance/` (schema
`qmsd.performance.v1`, `qmsd perf`), protocol assessments in `qmsd/data/performance/assess/`
(`qmsd assess --store`). `record_dict(record)["performance_record"]` returns the stored exact
record when one exists.

### `geometry` (optional)

Whatever the search family reasons about: `residual_checks` (dim A, the number of checks a
non-output column sees), `class_block_rank` and `plane_distance` for weighted supports,
`max_outputs_per_line` and `line_distance` for the m = 2 line-spread codes. These are the
quantities in `WEIGHTED_SUPPORTS.md` and `COVERAGE.md`; they are not recomputed by the tests.

## Status vocabulary (derived)

| status | certificate | examples |
|---|---|---|
| confirmed | exact (bounds meet) | the paper's search codes; [[1756,260,7]]_7; [[424,105,6]]_23; the RS family (closed-form MDS, re-checked by MITM for p <= 41) |
| partial | lower bound proved, exact value open; or closed form only | [[1746,270,>=7]]_7 (d = 7 or 8) |
| candidate | no rigorous bound | the old p = 17 / p = 19 d = 7 candidates |
| refuted | a lighter word was found later | the qutrit m = 7 cap code |

## Adding a record

Write a `qmsd.search-results.v1` bundle (the search CLI does this) and run
`qmsd catalog import bundle.json`. Each code may carry `construction`, `certificate`,
`performance`, `geometry`; missing blocks are inferred from the v1 fields (puncture columns ->
`punctured_rm`; `distance_certified` -> an exact engine certificate). A record cannot be imported
as `confirmed` unless its certificate is exact. For weighted supports the bundle must include the
`weighted_rm` construction with exactly k `outputs_in_X` (see
`runs/weighted/p7_m4_sphere_1756_260_bundle.json` and `runs/certify/p7_k270_seed407_bundle.json`).

`qmsd catalog show <artifact_id>` prints the record with its blocks; `--json` dumps everything,
including the stored performance record.

## Migration record

`qmsd catalog migrate` upgraded the 245 v1 files on 2026-10-29 (every v1 field kept; blocks inferred).
`scratchpad/catalog_enrich_v2.py` then filled the blocks whose history is known:

* the weighted p = 7 codes: `weighted_rm` constructions, staged certificates from
  `runs/certify/p7_k270_seed407.json` and the [[1756,260,7]]_7 proof, acceptance Monte Carlo;
* the Reed-Solomon family: closed-form certificates (exact; MITM re-check noted for p <= 41);
* the m = 2 line-spread codes at p = 19, 23: weight-6 line witnesses recomputed (exact MITM on the
  most-punctured lines; no weight <= 5 exists on the full code, so any weight-6 line word is the
  witness) and stored, geometry
  (`max_outputs_per_line`, `line_distance`), acceptance Monte Carlo.

Records imported earlier from the search engine keep the engine's exact MITM certificate.
