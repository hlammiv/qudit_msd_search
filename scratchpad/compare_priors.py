import json
from scipy.stats import chi2
rows=[r for r in map(json.loads, open("runs/surface/escape_keep_rules.jsonl")) if r.get("seed_off")==100000 and r["variant"]=="one+gap-tick1-closed"]
for r in rows:
    print("d",r["d"],"priors",bool(r.get("priors")))
    w=(r["accepted"]-r["nontrivial"])/r["n_clean"]
    for c in range(4,15):
        v=r["rules_syndrome_gap"][f"none:{c},{c}"]
        n=v["kept_res"]+w*v["kept_clean"]
        if not n: continue
        ev=v["fails_res"]+v["fails_clean"]
        f=v["fails_res"]+w*v["fails_clean"]
        kept=v["kept_res"]+v["kept_clean"]
        up=chi2.ppf(.9,2*ev+2)/2/kept
        print(f"  cut {c:2d} err {f/n:.2e} ev {v['fails_res']}+{v['fails_clean']} up90 {up:.2e} att/kept {r['shots']/n:7.1f}")
