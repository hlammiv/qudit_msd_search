# validation of branch_t_cult.py (2026-11-23); usage in MEAS_DISTILL 5k
import sys, numpy as np, time
import exact_t_cult as E
import branch_t_cult as BT
d=2; p=float(sys.argv[1]); N=int(sys.argv[2]); K=int(sys.argv[3])
cc, plan = BT.colorcult_plan(d, p, 5, "par-hooks-f1")
tl_q = {k:{plan.data[j]:t for j,t in v.items()} for k,v in plan.tl.items()}
t=time.time(); H = E.Hybrid(cc.c, cc.data, tl_q, BT.escaping_carrier(d)); print("dense build",round(time.time()-t),"s",flush=True)
n=cc.n; grid=H.grid
lab_z=((cc.MZ@grid)%d).T@(d**np.arange(cc.MZ.shape[0])); lab_x=((cc.H@grid)%d).T@(d**np.arange(cc.H.shape[0]))
emax=0
for a in range(d):
  for b in range(d):
    y=H.pauli_flat(H.out0[None,:],(a*plan.xl)[None,:]%d,(b*plan.zl)[None,:]%d); emax=max(emax,abs(np.vdot(H.out0,y[0])-plan.e[a,b]))
print("e dev",emax,flush=True)
rng=np.random.default_rng(9)
fp=BT.frame_pass(cc.c, plan.data, set(plan.flush), N, rng)
act=np.any(fp["flips"],axis=1)
for k,(xr,zr) in fp["rec"].items(): act|=np.any(xr,1)|np.any(zr,1)
ok=np.ones(N,bool)
for m in plan.trivial: ok&=fp["flips"][:,m]==0
ia=np.nonzero(act&ok)[0]
# prefer shots with X frames on data at T layers (branching)
W=np.zeros(N,int)
for k in plan.tl: W+=(fp["rec"][k][0]!=0).sum(1)
okp,_,_=BT.proxy_outcome(cc,plan,fp); ev0,_=plan.run_rows(fp, ia); accd=ia[ev0[:,0]>1e-9]; ia=np.concatenate([accd[W[accd]>=1][:K], accd[W[accd]==0][:K//4]])
print("selected W", W[ia][:20].tolist(), flush=True)
ev,mr=plan.run_rows(fp, ia)
res=[]
t=time.time()
for i in range(0,len(ia),4):
    sel=ia[i:i+4]
    rec={k:(xr[sel],zr[sel]) for k,(xr,zr) in fp["rec"].items()}
    psi=H.dense(rec, fp["flips"][sel], H.k0)
    w=(np.abs(psi)**2).sum(1)
    Ft=np.abs(psi@H.out0.conj())**2
    y0=psi*(lab_z==0)[None,:]
    y0f=np.fft.fftn(y0.reshape((len(sel),)+(d,)*n),axes=tuple(range(1,n+1)),norm="ortho").reshape(len(sel),-1)*(lab_x==0)[None,:]
    wt=(np.abs(y0f)**2).sum(1)
    res.append(np.stack([w,wt,Ft],1))
res=np.concatenate(res)
dev=np.max(np.abs(res-ev[:,[0,2,3]]))
print(f"D=5 d=2: {len(ia)} shots, max dev {dev:.2e}; w>0: {np.sum(res[:,0]>1e-12)}; sum w {res[:,0].sum():.4f} vs {ev[:,0].sum():.4f}; triv infid {np.sum(res[:,1]-res[:,2]):.4e} vs {np.sum(ev[:,2]-ev[:,3]):.4e}; max rows {mr}; dense {time.time()-t:.0f}s")
