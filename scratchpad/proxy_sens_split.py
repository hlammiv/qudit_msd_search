"""Cultivation-proxy sensitivity (MEAS_DISTILL 5k): split the closed-window escape error into residual-frame (cultivation)
and clean-frame (escape) failures, scale the residual part by X = 1, 3, 10, 30; prints p, cut, X, attempts/kept, central, UL
(UL = X UL(res) + w UL(clean)).  python scratchpad/proxy_sens_split.py > runs/surface/proxy_sens/scaled_p_in.txt"""
import json, sys
from scipy.stats import chi2
L=[json.loads(l) for l in open("runs/surface/escape_keep_rules.jsonl")]
ul=lambda n: chi2.ppf(0.9,2*(n+1))/2
sets={"1e-3":[12,15,17],"7e-4":[22],"5e-4":[20]}
for lab,ids in sets.items():
  for cut in range(5,12):
    k=f"none:{cut},{cut}"; fr=kr=fc=kc=0.0; nfc=0; shots=0; wmax=0
    for i in ids:
        r=L[i]; g=r["rules_syndrome_gap"][k]; w=(r["accepted"]-r["nontrivial"])/r["n_clean"]
        fr+=g["fails_res"]; kr+=g["kept_res"]; fc+=w*g["fails_clean"]; kc+=w*g["kept_clean"]; nfc+=g["fails_clean"]; shots+=r["shots"]
    den=kr+kc; att=shots/den
    wbar=fc/max(nfc,1e-9) if nfc else (sum(((L[i]["accepted"]-L[i]["nontrivial"])/L[i]["n_clean"]) for i in ids)/len(ids))
    for X in (1,3,10,30):
        cen=(X*fr+fc)/den; up=(X*ul(int(fr))+wbar*ul(nfc))/den
        print(lab,cut,X,f"{att:.2f}",f"{cen:.3e}",f"{up:.3e}")
