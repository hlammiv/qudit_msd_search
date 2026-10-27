"""Assign Z_{2^k} signs to the grafted / final codes extracted by graft_extract.py.  GF(2) system: every X/Z overlap of
size 2k has k sign disagreements (so the Z_4 overlap sum vanishes; the solution found vanishes exactly, hence all 2^k);
colour plaquettes keep the honeycomb tau; kept / extended / refined checks agree with their parents up to a unit; the
logicals Z^tau x Z^v (v on fresh |0>) and X x X^u (u on fresh |+>) commute with every grafted and final check.

python -m primepower_msd.graft_signs_solve <dir with graft_d5.json>
"""
import json, numpy as np, itertools as it, collections, math
import sys
S = sys.argv[1] if len(sys.argv) > 1 else '.'   # directory holding graft_d5.json; writes graft_d5_signed.json there
J=json.load(open(S+'/graft_d5.json'))
coords={int(k):tuple(v) for k,v in J['coords'].items()}
def stabs(key):
    out=[]
    for v in J[key]:
        if not v['data']: continue
        out.append((v['kind'][1:] or 'Z', frozenset(int(q) for q in v['data']), v['anc']))
    return out
C=stabs('colour'); G=stabs('grafted'); F=stabs('final')
cq=sorted(set(q for _,s,_ in C for q in s))
Zpl=[s for t,s,_ in C if t=='Z']
adj={q:set() for q in cq}
for a,b in it.combinations(cq,2):
    if math.dist(coords[a],coords[b])<1.5 and any(a in s and b in s for s in Zpl): adj[a].add(b); adj[b].add(a)
tau={cq[0]:1}; stack=[cq[0]]
while stack:
    i=stack.pop()
    for j in adj[i]:
        if j not in tau: tau[j]=-tau[i]; stack.append(j)
        assert tau[j]==-tau[i], 'not bipartite'
assert len(tau)==19, ('disconnected', len(tau))
assert all(sum(tau[q] for q in s)==0 for s in Zpl), 'tau not balanced'
print('tau sum', sum(tau.values()))
pl=Zpl; n=len(pl)
share=[[len(pl[i]&pl[j])>=2 for j in range(n)] for i in range(n)]
col=[-1]*n
def bt(i):
    if i==n: return True
    for c in range(3):
        if all(not(share[i][j] and col[j]==c) for j in range(i)):
            col[i]=c
            if bt(i+1): return True
    col[i]=-1; return False
assert bt(0); colours={pl[i]:col[i] for i in range(n)}
print('plaquette colours', sorted((sorted(coords[q] for q in s), c) for s,c in colours.items()))
inc={}
def var(key):
    if key not in inc: inc[key]=len(inc)
    return inc[key]
rows=[]; rhs=[]
def add_eq(vars_, b): rows.append(list(vars_)); rhs.append(b%2)
for code,L in (('G',G),('F',F)):
    for i,(t,s,_) in enumerate(L):
        for q in s: var((code,i,q))
    Xs=[(i,s) for i,(t,s,_) in enumerate(L) if t=='X']; Zs=[(i,s) for i,(t,s,_) in enumerate(L) if t=='Z']
    for i,sa in Xs:
        for j,sb in Zs:
            ov=sa&sb
            if not ov: continue
            assert len(ov)%2==0, ('odd overlap', code, len(ov))
            add_eq([inc[(code,i,q)] for q in ov]+[inc[(code,j,q)] for q in ov], len(ov)//2)

# ---- logical representatives (GF(2)): Z_L = Z on colour ∪ v (v ⊂ fresh |0>), commuting with all X checks of G and F;
#      X_L = X on colour ∪ u (u ⊂ fresh |+>), commuting with all Z checks of G and F.
fresh={int(q):b for q,b in J['fresh'].items()}
FZ=sorted(q for q,b in fresh.items() if b=='R'); FX=sorted(q for q,b in fresh.items() if b=='RX')
def gf2_solve(Arows, b):
    """solve A x = b over GF(2); Arows list of lists of column indices; returns one solution (min pivots) or None"""
    ncol=max((c for r in Arows for c in r), default=-1)+1
    M=np.zeros((len(Arows), ncol+1), dtype=np.uint8)
    for i,(r,bb) in enumerate(zip(Arows,b)):
        for c in r: M[i,c]^=1
        M[i,ncol]=bb%2
    rr=0; pv=[]
    for c in range(ncol):
        p=np.nonzero(M[rr:,c])[0]
        if len(p)==0: continue
        p=p[0]+rr; M[[rr,p]]=M[[p,rr]]
        msk=M[:,c].copy(); msk[rr]=0; M[msk==1]^=M[rr]; pv.append(c); rr+=1
        if rr==len(Arows): break
    z=M[rr:,:ncol].sum(axis=1)==0
    if np.any(M[rr:,ncol][z]==1): return None
    x=np.zeros(ncol,dtype=np.uint8)
    for k in range(rr-1,-1,-1):
        c=pv[k]; x[c]=(M[k,ncol]+int(M[k,c+1:ncol]@x[c+1:ncol]))%2
    return x
if 'LOGICALS' not in globals():
    Xall=[s for L in (G,F) for t,s,_ in L if t=='X']; Zall=[s for L in (G,F) for t,s,_ in L if t=='Z']
    colZ={q:i for i,q in enumerate(FZ)}; colX={q:i for i,q in enumerate(FX)}
    rowsv=[[colZ[q] for q in s if q in colZ] for s in Xall]; bv=[len([q for q in s if q in tau])%2 for s in Xall]
    v=gf2_solve(rowsv,bv); assert v is not None, 'no Z_L extension'
    rowsu=[[colX[q] for q in s if q in colX] for s in Zall]; bu=[len([q for q in s if q in tau])%2 for s in Zall]
    u=gf2_solve(rowsu,bu); assert u is not None, 'no X_L extension'
    ZL=set(tau)|{q for q in FZ if v[colZ[q]]}; XL=set(tau)|{q for q in FX if u[colX[q]]}
    print('Z_L support', len(ZL), 'fresh part', sorted(coords[q] for q in ZL-set(tau)))
    print('X_L support', len(XL), 'fresh part', sorted(coords[q] for q in XL-set(tau)))
    LOGICALS=True
# sign bits of the logical legs on fresh qubits are unknowns; colour legs fixed (Z: tau, X: +1)
def lbit(kind,q):
    if q in tau: return None, (int(tau[q]==-1) if kind=='Z' else 0)
    return var(('L',kind,q)), 0
for code,L in (('G',G),('F',F)):
    Xs=[(i,s) for i,(t,s,_) in enumerate(L) if t=='X']; Zs=[(i,s) for i,(t,s,_) in enumerate(L) if t=='Z']
    for i,sa in Xs:
        ov=[q for q in sa if q in ZL]
        if ov:
            assert len(ov)%2==0, ('odd Z_L overlap', code, sorted(coords[q] for q in ov))
            vs=[inc[(code,i,q)] for q in ov]; const=len(ov)//2
            for q in ov:
                vv,c=lbit('Z',q)
                if vv is None: const+=c
                else: vs.append(vv)
            add_eq(vs, const)
    for j,sb in Zs:
        ov=[q for q in sb if q in XL]
        if ov:
            assert len(ov)%2==0, ('odd X_L overlap', code, sorted(coords[q] for q in ov))
            vs=[inc[(code,j,q)] for q in ov]; const=len(ov)//2
            for q in ov:
                vv,c=lbit('X',q)
                if vv is None: const+=c
                else: vs.append(vv)
            add_eq(vs, const)
Cd={(t,s) for t,s,_ in C}
colour_bits={}
for t,s,_ in C:
    for q in s: colour_bits[(t,s,q)] = 0 if t=='X' else int(tau[q]==-1)
rel=0
for i,(t,s,_) in enumerate(G):
    sc=frozenset(q for q in s if q in tau)
    if sc and (t,sc) in Cd:
        g=var(('g','G',i))
        for q in sc: add_eq([inc[('G',i,q)], g], colour_bits[(t,sc,q)])
        rel+=1
    elif sc: print('  grafted stab with colour support not a colour stab:', t, sorted(coords[q] for q in s))
print('grafted stabs tied to colour stabs', rel)
Fset={(t,s):i for i,(t,s,_) in enumerate(F)}; Gset={(t,s):i for i,(t,s,_) in enumerate(G)}
same=0; prod=0
for i,(t,s,_) in enumerate(G):
    if (t,s) in Fset:
        j=Fset[(t,s)]; g=var(('g','GF',i))
        for q in s: add_eq([inc[('G',i,q)], inc[('F',j,q)], g], 0)
        same+=1
    else:
        parts=[(j,sf) for (tf,sf),j in Fset.items() if tf==t and sf<=s]
        found=False
        for k in range(2,4):
            for combo in it.combinations(parts,k):
                sets=[sf for _,sf in combo]
                if sum(len(x) for x in sets)==len(s) and frozenset().union(*sets)==s:
                    g=var(('g','GFp',i))
                    for j,sf in combo:
                        gj=var(('g','GFp',i,j))
                        for q in sf: add_eq([inc[('G',i,q)], inc[('F',j,q)], g, gj], 0)
                    prod+=1; found=True; break
            if found: break
        if not found: print('  grafted stab unaccounted in final:', t, sorted(coords[q] for q in s))
print('G stabs identical in F', same, 'G = product of F stabs', prod)
for i,(t,s,_) in enumerate(F):
    if (t,s) not in Gset: print('  new in final:', t, sorted(coords[q] for q in s))
nv=len(inc); m=len(rows)
A=np.zeros((m,nv+1),dtype=np.uint8)
for r,(vs,b) in enumerate(zip(rows,rhs)):
    for v in vs: A[r,v]^=1
    A[r,nv]=b
r=0; piv=[]
for ccol in range(nv):
    p=np.nonzero(A[r:,ccol])[0]
    if len(p)==0: continue
    p=p[0]+r; A[[r,p]]=A[[p,r]]
    mask=A[:,ccol].copy(); mask[r]=0
    A[mask==1]^=A[r]
    piv.append(ccol); r+=1
    if r==m: break
zero_rows=A[r:,:nv].sum(axis=1)==0
incons=bool(np.any(A[r:,nv][zero_rows]==1))
print('variables', nv, 'equations', m, 'rank', r, 'inconsistent', incons)
if not incons:
    x=np.zeros(nv,dtype=np.uint8)
    for k in range(r-1,-1,-1):
        ccol=piv[k]; x[ccol]=(A[k,nv]+int(A[k,ccol+1:nv]@x[ccol+1:nv]))%2
    def signed(code,L):
        return [dict(type=t, anc=a, legs={int(q):(-1 if x[inc[(code,i,q)]] else 1) for q in s}) for i,(t,s,a) in enumerate(L)]
    ZLs={int(q):(tau[q] if q in tau else (-1 if x[inc[('L','Z',q)]] else 1)) for q in ZL}
    XLs={int(q):(1 if q in tau else (-1 if x[inc[('L','X',q)]] else 1)) for q in XL}
    out=dict(coords={k:list(v) for k,v in coords.items()}, tau={int(q):v for q,v in tau.items()}, colours=[(sorted(s),c) for s,c in colours.items()], ZL=ZLs, XL=XLs,
             colour=[dict(type=t, anc=a, legs={int(q):(1 if t=='X' else tau[q]) for q in s}) for t,s,a in C],
             grafted=signed('G',G), final=signed('F',F), fresh=J['fresh'])
    json.dump(out, open(S+'/graft_d5_signed.json','w'))
    for name,L in (('colour',out['colour']),('grafted',out['grafted']),('final',out['final'])):
        bad=0; exact0=0; tot=0
        X=[v for v in L if v['type']=='X']; Z=[v for v in L if v['type']=='Z']
        for a in X:
            for b in Z:
                ov=set(a['legs'])&set(b['legs'])
                if ov:
                    ssum=sum(a['legs'][q]*b['legs'][q] for q in ov); tot+=1; bad+= ssum%4!=0; exact0+= ssum==0
        # logicals
        badL=sum(sum(a['legs'][q]*ZLs[q] for q in a['legs'] if q in ZLs)%4!=0 for a in X)+sum(sum(b['legs'][q]*XLs[q] for q in b['legs'] if q in XLs)%4!=0 for b in Z)
        print(name, 'overlapping pairs', tot, 'bad mod4', bad, 'exact zero', exact0, 'logical violations', badL)
