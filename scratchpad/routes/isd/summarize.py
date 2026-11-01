"""Print the numbers of the saved result JSONs (for NOTES.md)."""
import json, os, sys
R = os.path.dirname(os.path.abspath(__file__)) + "/runs"
def L(n): return json.load(open(os.path.join(R, n)))
if os.path.exists(R + "/p7a_k272_seed3_stern.json"):
    a = L("p7a_k272_seed3_stern.json")
    print("p7a", a["n"], a["checks"], a["E_law"])
    for r in a["runs"]: print("  run", r["found"], r["support"], r["iterations"], r["seconds"], r["params"]["ell"], r["params"]["z1"], r["p_success"]["6"])
    c = a["census_240s"]; print("  census iters", c["iterations"], "distinct", len(c["distinct_weight6_words"]), "hits", sum(c["hits_per_word"].values()), "P", c["p_success_6"])
if os.path.exists(R + "/p7b_k270_seed19_stern.json"):
    b = L("p7b_k270_seed19_stern.json")
    print("p7b", b["n"], b["checks"], b["E_law"])
    print("  w6", {k: b["w6"][k] for k in ("found", "iterations", "seconds", "p_success", "p_miss_one_word")})
    print("  w7 first", {k: b["w7_first"][k] for k in ("found", "support", "iterations", "seconds")})
    w = b["w7_census"]; print("  w7 census", {k: w[k] for k in ("iterations", "seconds", "total_hits", "distinct_words", "p_success_7", "N7_estimate_hits_over_itsP", "E7_law")})
if os.path.exists(R + "/p11_k2290_stern.json"):
    d = L("p11_k2290_stern.json"); print("p11", d["n"], d["checks"], d["E_law"])
    for k, r in d["runs"].items(): print("  w", k, {kk: r[kk] for kk in ("found", "iterations", "seconds", "p_success", "p_all_missed_model")}, r["params"])
    for k, rows in d["cost_table"].items():
        print("  cost table w", k); [print("    ", x) for x in rows]
for f in ("law_restricted_c13_s100.json", "law_restricted_c14_s120.json"):
    if os.path.exists(R + "/" + f):
        s = L(f); print(f, s.get("E_law"), s.get("summary"), [(x["rank"], x["lowest_weight"], x["n_words"], x["seconds"]) for x in s["subsets"]])
if os.path.exists(R + "/law_part2.json"):
    s = L("law_part2.json"); print("part2 summary", s.get("summary"))
    for n, r in s["candidates"].items():
        print("  ", n, r.get("checks"), r.get("E_law", {}).get("6"), r.get("P_no_w6_law"), r.get("stern6", {}).get("found"), r.get("stern6", {}).get("seconds"), r.get("stern6", {}).get("p_success", {}).get("6"), r.get("exact4"))
if os.path.exists(R + "/symmetry_check_k272.json"):
    for r in L("symmetry_check_k272.json"): print("sym", r["label"], r["n"], r["checks"], r["iterations"], r["hits"], r["distinct"], round(r["N6_estimate"], 1), round(r["E6_law"], 2), r["min_weight_seen"])
