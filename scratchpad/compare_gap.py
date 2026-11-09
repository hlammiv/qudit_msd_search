"""Syndrome-gap-field rows (closed or syndrome-only window) of escape_keep_rules.jsonl, pooled over replicas: error, 90% UL, attempts/kept.
    python -I scratchpad/compare_gap.py p variant [d]"""
import json, sys
from scipy.stats import chi2
p, var = float(sys.argv[1]), sys.argv[2]; d = int(sys.argv[3]) if len(sys.argv) > 3 else 4
sel = [r for r in map(json.loads, open("runs/surface/escape_keep_rules.jsonl")) if r["variant"] == var and r["p"] == p and (r.get("d") or 4) == d]
print(var, p, d, len(sel), "runs", sum(r["shots"] for r in sel), "shots")
for c in (0, 4, 6, 7, 8, 9, 10, 12, 14):
    n = f = er = ec = kept = shots = 0
    for r in sel:
        v = r["rules_syndrome_gap"][f"none:{c},{c}"]; w = (r["accepted"] - r["nontrivial"]) / r["n_clean"]
        n += v["kept_res"] + w * v["kept_clean"]; f += v["fails_res"] + w * v["fails_clean"]
        er += v["fails_res"]; ec += v["fails_clean"]; kept += v["kept_res"] + v["kept_clean"]; shots += r["shots"]
    if not kept:
        continue
    up = chi2.ppf(.9, 2 * (er + ec) + 2) / 2 * (n / kept) / n
    print(f"  cut {c:2d} err {f/n:.2e} ev {er}+{ec} up90 {up:.2e} att/kept {shots/n:7.1f}")
