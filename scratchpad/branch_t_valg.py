# validation of branch_t_cult.py (2026-11-23); usage in MEAS_DISTILL 5k
import sys, numpy as np, time
import branch_t_cult as BT
d=int(sys.argv[1]); p=float(sys.argv[2]); N=int(sys.argv[3]); proxy=sys.argv[4]=="proxy"
t=time.time(); g, plan = BT.growcult_plan(d, p, "par-hooks-f1", proxy=proxy); print("build",round(time.time()-t,1),"s; flush",len(plan.flush),"random",len(plan.rnd), flush=True)
print("e table:", np.round(plan.e,4).tolist())
# noiseless
fp0 = BT.frame_pass(g.c, plan.data, set(plan.flush), 4, np.random.default_rng(0), scale=0.0)
ev0,_ = plan.run_rows(fp0, np.arange(4)); print("noiseless", ev0.tolist(), flush=True)
rng=np.random.default_rng(11)
fp = BT.frame_pass(g.c, plan.data, set(plan.flush), N, rng)
okp, fail, nontriv = BT.proxy_outcome(g, plan, fp)
t=time.time(); ev, mr = plan.run_rows(fp, np.arange(N)); dt=time.time()-t
w=ev[:,0]; F=ev[:,1]
print(f"proxy acc {okp.sum()} fails {fail.sum()} triv-fails {(fail&~nontriv).sum()}; engine sum w {w.sum():.4f} infid {np.sum(w-F):.4f} triv w {ev[:,2].sum():.3f} triv infid {np.sum(ev[:,2]-ev[:,3]):.4f}; max rows {mr}; {dt:.1f}s for {N}")
if proxy:
    pf=np.zeros(N); pf[np.nonzero(okp)[0][fail]]=1
    print("shot-by-shot: accept mismatch", np.sum(np.abs(w-okp)>1e-9), " fail mismatch", np.sum(np.abs((w-F)-pf)>1e-9))
    bad=np.nonzero(np.abs(w-okp)>1e-9)[0][:5]; print("examples", [(int(i), float(w[i]), bool(okp[i])) for i in bad])
