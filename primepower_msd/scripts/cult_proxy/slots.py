import sys, collections
d=int(sys.argv[1]); NS=int(sys.argv[2]); sys.argv=['x',str(d)]
exec(open('proxy_validate.py').read().split('out = {')[0])
sups=set()
for sup in it.combinations(range(n),3):
    for Hc in (MZ,H):
        for v in it.product(range(1,d),repeat=3):
            e=np.zeros(n,dtype=np.int64); e[list(sup)]=v
            if not ((Hc@e)%d).any(): sups.add(sup)
sups=sorted(sups); print('logical supports',len(sups))
rng=np.random.default_rng(0)
tot=collections.defaultdict(lambda:[0.,0.,0.,0.])
pats=list(it.product(range(4),repeat=3))
for sup in [sups[j] for j in rng.choice(len(sups),min(NS,len(sups)),replace=False)]:
    for sl in pats:
        key=tuple(sorted(sl))
        for pv in it.product(paulis,repeat=3):
            fs=[(s,i,a,b) for s,i,(a,b) in zip(sl,sup,pv)]
            aT,uT=run("T",fault_set(fs)); aS,uS=run("S",fault_set(fs))
            t=tot[key]; t[0]+=uT; t[1]+=uS; t[2]+=aT; t[3]+=aS
for k,v in sorted(tot.items()): print(k, 'u_T %.3f u_S %.3f  acc_T %.2f acc_S %.2f'%tuple(v))
print('TOTAL u_T %.4f u_S %.4f'%(sum(v[0] for v in tot.values()), sum(v[1] for v in tot.values())))
