"""Pool syndrome-gap rows of escape_keep_rules.jsonl over replicas (same p, d, variant) and print p_total with a 90% upper limit."""
import json, sys
from scipy.stats import chi2
p, d, var = float(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
rs = [json.loads(l) for l in open(sys.argv[4] if len(sys.argv) > 4 else "runs/surface/escape_keep_rules.jsonl")]
sel = [r for r in rs if r.get("variant") == var and r["p"] == p and r.get("d", 4) == d]
print(len(sel), "replicas", sum(r["shots"] for r in sel), "shots")
for rule in ("rules_syndrome_gap", "rules_readout_fine"):
    print(rule)
    for c in range(0, 17):
        k = f"none:{c},{c}"
        n = fails = 0.0; kc = kr = fr = fc = 0; shots = 0
        for r in sel:
            v = r[rule][k]; w = (r["accepted"] - r["nontrivial"]) / r["n_clean"]
            n += v["kept_res"] + w * v["kept_clean"]; fails += v["fails_res"] + w * v["fails_clean"]
            kr += v["kept_res"]; fr += v["fails_res"]; kc += v["kept_clean"]; fc += v["fails_clean"]; shots += r["shots"]
        wbar = fails / max(fr + fc, 1) if fr + fc else n / max(kr + kc, 1)
        ev = fr + fc
        up = chi2.ppf(0.9, 2 * ev + 2) / 2 * (n / max(kr + kc, 1)) / n   # 90% upper on events, ~ weight per event
        print(f"  cut {c:2d} kept_res {kr:6d} fails_res {fr:3d} kept_clean {kc:7d} fails_clean {fc:3d} p_total {fails/n:.3e} up90~{up:.2e} att/kept {shots/n:8.1f}")
