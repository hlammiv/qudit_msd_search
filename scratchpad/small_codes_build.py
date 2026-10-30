"""Small effective codes for every prime <= 97: (1) import the whole Reed-Solomon family into the
catalog (construction reed_solomon, closed-form exact certificate) with exact performance records
(engine: mds closed form); (2) write SMALL_CODES.md: per prime, the RS codes with their exact
threshold / leading coefficient / acceptance, the m = 2 search codes in the catalog, and the
literature k = 1 codes at p = 3, 5, 7."""
import json, math, collections, glob, datetime
from pathlib import Path
from qmsd.catalog import construction_block, certificate_block, performance_block, import_search_export, catalog_directory
from qmsd.performance import write_performance, PERF_SCHEMA
from qmsd.rs_family import r_max_rs
from qmsd.results import load_result_catalog
rows = json.load(open("runs/coverage/rs_family_all.json"))["rows"]
today = datetime.date.today().isoformat()
# ---- 1. catalog + performance records for all 329 RS codes
codes = []
for r in rows:
    p, n, k, d, rr = r["p"], r["n"], r["k"], r["d"], r["r"]
    cert = certificate_block(d, d, lower_method="closed-form MDS distance d = r - k + 2 of the punctured Reed-Solomon code (RESULTS.md sec 4; shortened/punctured MDS codes are MDS)",
                             upper_method="closed-form MDS: any d columns of the [n, r+1-k] shortened RS dual are ... explicit: the MDS code has words of every weight >= d", machine="closed form", date=today)
    cert["exact"] = True
    codes.append({"p": p, "m": 1, "r": rr, "n": n, "k": k, "d": d, "gamma": r["gamma"], "A_d": None, "A_d_exact": None, "full_rank": True, "distance_certified": True,
                  "label": f"[[{n},{k},{d}]]_{p}", "puncture_columns_1indexed": list(range(1, k + 1)), "w": None,
                  "construction": construction_block("reed_solomon", p, 1, rr, puncture_columns_1indexed=list(range(1, k + 1))),
                  "certificate": cert, "performance": performance_block(exact=True, notes="exact closed-form MDS performance (qmsd.rs_family); record in qmsd/data/performance"),
                  "note": f"Reed-Solomon (m = 1) family member, k = {k} of r = {rr}; puncture set is immaterial (puncture-invariant). Exact threshold {r['threshold_avg']:.4f} (per output), c = {r['c_eps_avg']:.3g}."})
Path("runs/coverage").mkdir(exist_ok=True)
json.dump({"schema": "qmsd.search-results.v1", "search": {"p": "5-97", "m": 1, "sampler": "reed_solomon_family_all_k (scratchpad/rs_family_table.py)", "trials": 0, "seed": 0, "target_k": None, "top": len(codes), "jobs": 1, "max_triples": 0, "min_keep_distance": None}, "codes": codes},
          open("runs/coverage/rs_family_all_bundle.json", "w"), indent=1)
res = import_search_export("runs/coverage/rs_family_all_bundle.json"); print("catalog import:", len(res["imported"]), "new,", len(res["unchanged"]), "already present")
# performance records (exact) for every RS record in the catalog
recs = {(r.p, r.n, r.k): r for r in load_result_catalog() if r.construction and r.construction.get("family") == "reed_solomon"}
byrow = {(r["p"], r["n"], r["k"]): r for r in rows}
nperf = 0
for key, rec in recs.items():
    r = byrow.get(key)
    if r is None: continue
    out = {"artifact_id": rec.artifact_id, "label": rec.label, "p": rec.p, "n": rec.n, "k": rec.k, "convention": "eps = paper (2610.08012); delta = p/(p-1) eps",
           "dim_A": r["dim_A"], "tier": 1, "engine_avg": "mds closed form (qmsd.rs_family)", "engine_block": "mds closed form (qmsd.rs_family)",
           "avg": {"d": r["d"], "L_d": r["L_d_avg"], "c_eps": r["c_eps_avg"], "c_delta": r["c_eps_avg"] * ((rec.p - 1) / rec.p) ** r["d"] * (rec.p - 1) / rec.p if False else None,
                   "threshold_eps": r["threshold_avg"], "threshold_delta": None if r["threshold_avg"] is None else r["threshold_avg"] * rec.p / (rec.p - 1), "threshold_status": "found" if r["threshold_avg"] else "none below (p-1)/p", "threshold_rigorous": True},
           "block": {"d": r["d"], "c_eps": r["c_eps_block"], "threshold_eps": r["threshold_block"], "threshold_delta": None if r["threshold_block"] is None else r["threshold_block"] * rec.p / (rec.p - 1), "threshold_status": "found" if r["threshold_block"] else "none below (p-1)/p", "threshold_rigorous": True},
           "P_acc_eps_0.01": r["P_acc_1e-2"], "P_acc_eps_0.001": r["P_acc_1e-3"], "status": "exact", "seconds": 0.0}
    write_performance(out); nperf += 1
print("performance records written:", nperf)
# ---- 2. SMALL_CODES.md
cat = list(load_result_catalog())
byp = collections.defaultdict(list)
for r in rows: byp[r["p"]].append(r)
lit = {3: [("[[15,1,3]]_3", "doubling (2610.29457)", "", ""), ("[[19,1,3]]_3", "cyclic B (2610.08012)", "5.25", "0.226"), ("[[8,1,2]]_3", "RM", "2.0", "0.211"), ("[[37,1,5]]_3", "doubling", "", ""), ("[[72,9,3]]_3", "paper Table 3", "", ""), ("[[206,37,4]]_3", "paper Table 3", "15.2 (per output)", "0.0616")],
       5: [("[[4,1,2]]_5", "RS", "1.5", "0.363"), ("[[11,1,3]]_5", "cyclic B", "1.44", "0.346"), ("[[16,1,4]]_5", "cyclic", "3.38", "0.363"), ("[[33,1,5]]_5", "doubling", "", ""), ("[[20,5,2]]_5", "paper Table 3", "", ""), ("[[112,13,3]]_5", "paper Table 3", "", "")],
       7: [("[[6,1,2]]_7", "RS", "2.5", "0.232"), ("[[7,1,3]]_7", "cyclic B", "0.97", "0.446"), ("[[12,1,3]]_7", "cyclic A", "1.11", "0.366"), ("[[20,1,4]]_7", "doubling", "", ""), ("[[44,5,3]]_7", "m = 2 search (catalog)", "1.11", "0.173")]}
L = ["# Small effective codes for every prime p <= 97", "",
     f"*{today}. Companion to `COVERAGE.md` (which ranks by gamma). Here the question is the practical one: at each prime,",
     "which small codes exist, how good must the raw magic states be (threshold eps*), how fast does one round suppress",
     "error (eps_out ~ c eps^d), and how many blocks are accepted (P_acc). All thresholds and coefficients are EXACT:",
     "Reed-Solomon codes from the closed-form MDS engine (`qmsd.rs_family`, `runs/coverage/rs_family_all.json`), the",
     "rest from the enumerating engine (`qmsd/data/performance/`). eps is the convention of arXiv:2610.08012 (probability",
     "of a nontrivial Z error per input).*", "",
     "## Reading guide", "",
     "- **k = 1 Reed-Solomon codes [[p-1, 1, r+1]]_p** are the best single-output codes at every prime >= 11: thresholds",
     "  0.5-0.7 and leading coefficients that fall like p^{-(d-1)}; e.g. [[22,1,8]]_23 has eps* = 0.64 and c = 1.3e-4.",
     "- Larger k trades distance for rate along the family (d = r - k + 2); the minimum-gamma member is the last column.",
     "- Below p = 11 the RS family is trivial (d = 2) and the small codes are the cyclic / doubling / paper codes.", "",
     "## Reed-Solomon family, exact (per-output threshold eps*, coefficient c of eps_out ~ c eps^d, P_acc at eps = 1e-2)", ""]
L += ["| p | r | k = 1 code | eps* | c | P_acc(1e-2) | min-gamma member | gamma | eps* | all k |", "|---|---|---|---|---|---|---|---|---|---|"]
for p in sorted(byp):
    rs = byp[p]; k1 = rs[0]; bg = min((r for r in rs if r["gamma"]), key=lambda r: r["gamma"]) if any(r["gamma"] for r in rs) else None
    allk = ", ".join(f"[[{r['n']},{r['k']},{r['d']}]] eps*={r['threshold_avg']:.2f}" for r in rs)
    L.append(f"| {p} | {k1['r']} | [[{k1['n']},1,{k1['d']}]] | {k1['threshold_avg']:.3f} | {k1['c_eps_avg']:.2e} | {k1['P_acc_1e-2']:.2f} | " +
             (f"[[{bg['n']},{bg['k']},{bg['d']}]] | {bg['gamma']:.3f} | {bg['threshold_avg']:.3f} |" if bg else "- | - | - |") + f" {allk} |")
L += ["", "## Search and literature codes with small blocks (p <= 29)", "",
      "| p | code | source | c (block unless noted) | eps* | notes |", "|---|---|---|---|---|---|"]
for p in (3, 5, 7):
    for lab, src, c, thr in lit[p]:
        L.append(f"| {p} | {lab} | {src} | {c or '-'} | {thr or '-'} | |")
m2 = [r for r in cat if r.construction and r.construction.get("family") == "punctured_rm" and r.m == 2 and r.status == "confirmed" and r.p >= 11]
best = {}
for r in m2:
    key = (r.p, r.k)
    if key not in best or (r.d or 0) > (best[key].d or 0): best[key] = r
for (p, k), r in sorted(best.items()):
    L.append(f"| {p} | {r.label} | m = 2 search (catalog `{r.artifact_id}`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |")
L += ["", "## How to extend", "",
      "- `scratchpad/rs_family_table.py` regenerates the RS table; `scratchpad/small_codes_build.py` rebuilds this file and the catalog/performance records.",
      "- m = 2 searches at small k (high d, n ~ 100-300) for p = 11, 13 run on loranne (`runs/coverage/small/`); p >= 17 at d >= 6 needs the 3+3 kernel memory (lenore).",
      "- Cost per output (raw states per accepted output at a target error) comes from `qmsd protocol` / `qmsd assess`; not yet tabulated here."]
open("SMALL_CODES.md", "w").write("\n".join(L) + "\n"); print("SMALL_CODES.md written")
