"""One-off enrichment after the v1 -> v2 migration (2026-10-29): fill the structured blocks of the
records whose history is known: the weighted p = 7 codes (construction + staged certificates),
the Reed-Solomon family (closed-form certificates), the m = 2 line-spread codes (recomputed
weight-6 line witnesses, geometry), and the acceptance Monte Carlo pointers."""
import json, glob, math, itertools, numpy as np
from pathlib import Path
from qmsd.catalog import construction_block, certificate_block, performance_block, status_from_certificate, import_search_export, catalog_directory
from qmsd.reedmuller import points, rm_generator
from qmsd.structured_ad import _rank, _null_space, mitm_streamed
from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.performance import canonical_logical_rows
cat = catalog_directory()
def load(aid): return json.load(open(cat / f"{aid}.json"))
def save(aid, rec): (cat / f"{aid}.json").write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n")
pacc = json.load(open("runs/certify/pacc_mc.json"))
def pacc_block(label):
    v = pacc.get(label); return None if v is None else {"source": "runs/certify/pacc_mc.json", "P_acc": {e: x["mc"] for e, x in v["P_acc"].items()}, "note": "Monte Carlo, 2-4e5 shots; equals (1-eps)^n within 1%"}
# ---- 1. weighted p = 7 records
rec407 = json.load(open("runs/certify/p7_k270_seed407.json"))
cons = construction_block("weighted_rm", 7, 4, 7, charge={"kind": "quadric", "diag": [1, 1, 1, 1], "lin": [0, 0, 0, 0], "const": 0}, outputs_in_X=rec407["cols_in_X"])
st = rec407["stages"]
stages = [{"stage": "block6", "result": "no class-block word of weight <= 6 (serial MITM <= 5, 3+3 scan for 6)", "seconds": st["block6"]["seconds"]},
          {"stage": "full5", "result": "no word of weight <= 5 on the full code (serial MITM <= 4, parallel 2+3 weight-5 scan)", "seconds": st["full5"]["seconds"]},
          {"stage": "mixed6", "result": f"no weight-6 word touching a non-class column ({st['mixed6']['line_word_candidates']} line-word candidates, 0 extensions)", "seconds": st["mixed6"]["seconds"]},
          {"stage": "witness7", "result": "plane distance 8; no planar or 3-flat weight-7 word; full 3+4 block scan deferred", "seconds": None}]
cert = certificate_block(7, None, lower_method="block6 + full5 + mixed6 (scratchpad/certify_candidate.py)", upper_method="open: d = 7 or 8 (no weight-7 witness found in planar / 3-flat scans)", stages=stages, machine="lenore, 8 threads nice 19", date="2026-10-29")
for f in glob.glob(str(cat / "imported-p7-n1746-k270-*.json")):
    aid = Path(f).stem; r = load(aid)
    r.update({"construction": cons, "certificate": cert, "status": status_from_certificate(cert), "family": "weighted support (sphere charge), sweep candidate",
              "performance": performance_block(acceptance_mc=pacc_block("[[1746,270,>=7]]_7"), notes="exact enumeration infeasible (dim A = 56)"),
              "geometry": {"residual_checks": rec407["G0_rows"], "class_block_rank": rec407["block_rank"], "class_cols": rec407["class_cols"], "plane_distance": rec407["d_plane"], "sweep": "runs/weighted/plane_sweep (k=270 seed 407)"}})
    save(aid, r); print("enriched", aid)
# the [[1756,260,7]]_7 code: new record with exact certificate + witness
r260 = json.load(open("runs/weighted/p7_m4_sphere_1756_260_cols.json"))
wit = {"columns": r260["weight7_witness_G0_cols"], "coefficients": r260["witness7_coeffs"][0], "weight": 7}
cert260 = certificate_block(7, 7, lower_method="no weight <= 5 (numba MITM + parallel weight-5 scan); no weight-6 block word (3+3); no mixed weight-6 word (RM(20,4) line lemma, scratchpad/mixed_w6_check.py)", upper_method="explicit weight-7 word (3+4 block scan)", witness=wit,
                            stages=[{"stage": "certificate", "result": "three-part proof + witness, see WEIGHTED_SUPPORTS.md", "seconds": None}], machine="lenore (2026-10-27)", date="2026-10-27")
code260 = {"p": 7, "m": 4, "r": 7, "n": 1756, "k": 260, "d": 7, "gamma": math.log(1756 / 260) / math.log(7), "A_d": None, "A_d_exact": None, "full_rank": True, "distance_certified": True,
           "label": "[[1756,260,7]]_7", "puncture_columns_1indexed": None, "w": None,
           "construction": construction_block("weighted_rm", 7, 4, 7, charge={"kind": "quadric", "diag": [1, 1, 1, 1], "lin": [0, 0, 0, 0], "const": 0}, outputs_in_X=r260["cols_in_X"]),
           "certificate": cert260, "performance": performance_block(acceptance_mc=pacc_block("[[1756,260,7]]_7"), notes="exact enumeration infeasible (dim A = 66); weight-7 census queued (CERTIFY_P7.md)"),
           "geometry": {"residual_checks": 66, "class_block_rank": 31, "plane_distance": 7},
           "note": "First gamma < 1 code at p = 7 (WEIGHTED_SUPPORTS.md): weighted sphere-charge support on F_7^4, outputs on one cube class; wm4_bracket.py seed 7, k = 260 sample 1.",
           "distance_evidence": "d = 7 exact: three-part lower bound + explicit weight-7 word (WEIGHTED_SUPPORTS.md)"}
json.dump({"schema": "qmsd.search-results.v1", "search": {"p": 7, "m": 4, "sampler": "weighted sphere-charge bracket (wm4_bracket.py)", "seed": 7, "target_k": 260, "trials": 2, "top": 1, "jobs": 1, "max_triples": 0, "min_keep_distance": None}, "codes": [code260]}, open("runs/weighted/p7_m4_sphere_1756_260_bundle.json", "w"), indent=1)
print(import_search_export("runs/weighted/p7_m4_sphere_1756_260_bundle.json"))
# ---- 2. Reed-Solomon family
for f in glob.glob(str(cat / "imported-p*-n*-k*.json")):
    aid = Path(f).stem; r = load(aid)
    if r.get("m") != 1: continue
    mitm = "MITM exact d = r - k + 2 (min_dependent_columns)" if r["p"] <= 41 else ""
    cert = certificate_block(r["d"], r["d"], lower_method="closed-form MDS distance d = r - k + 2 of the punctured Reed-Solomon code (RESULTS.md sec 4)" + (f"; {mitm}" if mitm else ""), upper_method="closed-form MDS (any k+... columns) " + ("+ MITM dependent set" if mitm else ""), machine="local" if mitm else "", date="2026-10-28")
    cert["exact"] = True
    r.update({"certificate": cert, "status": "confirmed", "construction": construction_block("reed_solomon", r["p"], 1, r["r"], puncture_columns_1indexed=r["puncture_columns"]),
              "performance": performance_block(exact=r["p"] <= 41, notes="exact curve + threshold stored for p <= 41 (qmsd perf)" if r["p"] <= 41 else "dim A too large for exact enumeration; closed-form A_d = C(n,d)(p-1)"),
              "distance_evidence": cert["lower_method"]})
    save(aid, r); print("RS", aid, r["status"])
# ---- 3. m = 2 line-spread codes (p = 19, 23): recompute the weight-6 line witness, store geometry
P2 = {}
for f in glob.glob(str(cat / "imported-p*-n*-k*.json")):
    aid = Path(f).stem; r = load(aid)
    if r.get("m") != 2 or "line_spread_sweep" not in r.get("note", ""): continue
    p, k, rr = r["p"], r["k"], r["r"]; cols = r["puncture_columns"]
    if p not in P2:
        P2[p] = (np.asarray(points(2, p), dtype=np.int64), np.asarray(rm_generator(rr, 2, p), dtype=np.int64) % p)
    P, G = P2[p]
    built = build_triorthogonal_code(p, 2, rr, cols, G=G); G0 = np.asarray(built["G0"], dtype=np.int64) % p
    cs = set(c - 1 for c in cols); keep = [l for l in range(len(P)) if l not in cs]; pos = {l: j for j, l in enumerate(keep)}
    lines = []
    for a in range(p):
        key = (P[:, 0] + a * P[:, 1]) % p; lines += [np.flatnonzero(key == c) for c in range(p)]
    lines += [np.flatnonzero(P[:, 1] == c) for c in range(p)]
    # no weight <= 5 word exists on the full code (certified), so any weight-6 line word is THE witness;
    # the lightest line words sit on the most-punctured lines: search those first, stop at the first hit
    jcount = [sum(1 for i in ln if i in cs) for ln in lines]; maxline = max(jcount)
    best = None
    for ln, j in sorted(zip(lines, jcount), key=lambda x: -x[1]):
        cc = np.array([pos[i] for i in ln if i in pos], dtype=np.int64)
        if len(cc) < 6: continue
        M = G0[:, cc]; M = M[[tt for tt in range(M.shape[0]) if M[tt].any()]]
        if _rank(M.T.copy(), p) == len(cc): continue
        ws = list(itertools.islice(mitm_streamed(M, 6, p), 1))
        if ws:
            supp, vals = ws[0]; best = (6, [int(cc[s]) for s in supp], [int(v) for v in vals]); break
    if best is None:
        print("NO weight-6 line word for", aid, "-> lower bound only (d >= 6)"); continue
    vec = np.zeros(G0.shape[1], dtype=np.int64); vec[best[1]] = best[2]; assert ((G0 @ vec) % p == 0).all()
    wit = {"columns": best[1], "coefficients": best[2], "weight": best[0]}
    cert = certificate_block(r["d"], r["d"], lower_method="no weight <= 5: serial numba MITM d <= 4 + thread-parallel 2+3 weight-5 scan (scratchpad/m2_line_sweep.py)", upper_method="explicit line-supported weight-6 word (exact MITM on the line)", witness=wit, machine="loranne / local / lenore", date="2026-10-28")
    r.update({"certificate": cert, "status": "confirmed", "construction": construction_block("punctured_rm", p, 2, rr, puncture_columns_1indexed=cols),
              "performance": performance_block(acceptance_mc=pacc_block(r["label"]), notes="exact enumeration infeasible (dim A >= 13)"),
              "geometry": {"residual_checks": int(G0.shape[0]), "max_outputs_per_line": int(maxline), "line_distance": best[0]},
              "distance_evidence": cert["lower_method"] + "; " + cert["upper_method"]})
    save(aid, r); print("line-sweep", aid, "witness weight", best[0], "max_line", maxline)
print("done")
