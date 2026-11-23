# validation of branch_t_cult.py (2026-11-23); usage in MEAS_DISTILL 5k
import sys, numpy as np, time
import exact_t_cult as E
import branch_t_cult as BT
d=int(sys.argv[1]); p=float(sys.argv[2]); N=int(sys.argv[3])
cc, plan = BT.colorcult_plan(d, p, 3, "par-hooks-f1")
tl = plan.tl
tl_q = {k:{plan.data[j]:t for j,t in v.items()} for k,v in tl.items()}
H = E.Hybrid(cc.c, cc.data, tl_q, BT.escaping_carrier(d))
dec = E.Decoder(cc, H.out0)
# e table check vs dense out0
out0 = H.out0.reshape((d,)*cc.n)
Dn = E.Dense(cc.n, d)
emax=0
for a in range(d):
  for b in range(d):
    y = H.out0[None,:]
    y = H.pauli_flat(y, (a*plan.xl)[None,:]%d, (b*plan.zl)[None,:]%d)
    ev = np.vdot(H.out0, y[0])
    emax=max(emax, abs(ev-plan.e[a,b]))
print("e table max dev", emax, flush=True)
rng=np.random.default_rng(5)
fp = BT.frame_pass(cc.c, plan.data, set(plan.flush), N, rng)
act = np.any(fp["flips"],axis=1)
for k,(xr,zr) in fp["rec"].items(): act |= np.any(xr,1)|np.any(zr,1)
ia=np.nonzero(act)[0][:int(sys.argv[4]) if len(sys.argv)>4 else 400]
t=time.time(); ev, mr = plan.run_rows(fp, ia); tb=time.time()-t
rec={k:(xr[ia],zr[ia]) for k,(xr,zr) in fp["rec"].items()}
t=time.time(); psi=H.dense(rec, fp["flips"][ia], H.k0); w=(np.abs(psi)**2).sum(1); Fr,wt,Ft=dec.fidelity(psi); td=time.time()-t
dev=np.max(np.abs(np.stack([w,Fr,wt,Ft],1)-ev))
print(f"d={d} D=3: {len(ia)} shots, max dev {dev:.2e}; nonzero w {np.sum(w>1e-12)}, sum w {w.sum():.3f} vs {ev[:,0].sum():.3f}; infid {np.sum(w-Fr):.4e} vs {np.sum(ev[:,0]-ev[:,1]):.4e}; max rows {mr}; branch {tb:.2f}s dense {td:.2f}s")
