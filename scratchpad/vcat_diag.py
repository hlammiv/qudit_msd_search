import sys, numpy as np
from primepower_msd.cult_circuit import ColorCult
from primepower_msd.zdsim import add_tick_idle
d=int(sys.argv[1]); p=float(sys.argv[2]); shots=int(sys.argv[3])
base=dict(rounds=1, checks_per_round=2, hooks=True, parallel=True, final_block=True)
for vc in (False, True):
  for f in ((1.0,) if len(sys.argv)>4 else (0.0,1.0)):
    cc=ColorCult(3,d,p,verify_cat=vc,**base)
    nloc=sum(1 for op in cc.c.ops if op[0] in("DEP1","DEP2","XERR"))
    if f: cc.c=add_tick_idle(cc.c,f*p,no_anchor=(cc.vanc or []) if len(sys.argv)>4 else ())
    r=cc.simulate(shots,seed=7,batch=100000)
    at=r["accepted_trivial_final_syndrome"]
    print(f"vcat={vc} f={f}: acc {r['accept_rate']:.4f} pL {r['p_L']:.3e} triv {r['failures_trivial_final_syndrome']/max(at,1):.3e} ({r['failures_trivial_final_syndrome']} ev)",flush=True)
