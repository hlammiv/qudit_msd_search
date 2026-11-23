"""Footprints for each line of proxy_sens_split output (cuts 5-11; fewest qudits, then fewest attempts).
python scratchpad/proxy_sens_footprint.py runs/surface/proxy_sens/scaled_p_in.txt"""
import subprocess, re, sys
rows=[l.split() for l in open(sys.argv[1])]
best={}
for lab,cut,X,att,cen,up in rows:
    pv=float(lab)
    for kind,q in (("UL",up),("cen",cen)):
        if float(q)==0: continue
        o=subprocess.run(["python3","scratchpad/footprint_tick_proxy.py","1",f"e2e={pv},{att},{q}","escape-gap","census-sparse"],capture_output=True,text=True,env={**__import__("os").environ,"PYTHONPATH":".:scratchpad"}).stdout
        tag=f"{pv:.0e}".replace("e-0","e-0")
        for tgt in ("1e-10","1e-12"):
            m=[l for l in o.splitlines() if l.startswith(f"p={pv:.0e}") and f"target={tgt}" in l]
            if not m: continue
            mm=re.search(r"fewest qudits D=(\d+) (.*?): (\d+) physical.*?([\d.]+) cultivation attempts/output",m[0])
            D,code,qd,ao=int(mm.group(1)),mm.group(2),int(mm.group(3)),float(mm.group(4))
            rounds=code.count("[[")
            key=(lab,X,kind,tgt)
            cand=(qd,ao,int(cut),D,rounds)
            if key not in best or cand<best[key]: best[key]=cand
for k in sorted(best, key=lambda k:(k[0],k[3],k[2],int(k[1]))):
    qd,ao,cut,D,r=best[k]; print(*k, f"cut {cut} D={D} rounds {r} qudits {qd} att/out {ao:.0f}")
