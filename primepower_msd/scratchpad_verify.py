import itertools, numpy as np
# Concrete [[5,2,3]]_4 non-CSS Z_4 code from settled context
A=[[1,2,1,1,2],[1,1,1,2,1],[3,2,0,1,2],[2,1,2,0,2]]
B=[[1,0,3,3,1],[2,1,3,0,2],[0,1,1,1,0],[2,1,3,1,3]]
A=np.array(A)%4; B=np.array(B)%4
n=5
def symp(a1,b1,a2,b2):
    return int((a1@b2 - a2@b1))%4
print("pairwise symplectic (need 0 mod4):")
ok=True
for i in range(4):
    for j in range(i+1,4):
        s=symp(A[i],B[i],A[j],B[j])
        print(i,j,s)
        if s!=0: ok=False
print("all commute:",ok)

# self-symplectic (order check) 
for i in range(4):
    print("self",i,symp(A[i],B[i],A[i],B[i]))

# Build stabilizer group from integer combos of the 4 generators mod4
gens=[(A[i],B[i]) for i in range(4)]
S=set()
for c in itertools.product(range(4),repeat=4):
    a=np.zeros(n,int); b=np.zeros(n,int)
    for k in range(4):
        a=(a+c[k]*gens[k][0])%4
        b=(b+c[k]*gens[k][1])%4
    S.add((tuple(a),tuple(b)))
print("|S| =",len(S))

# symplectic weight of a Pauli (a|b): number of coords where (a_i,b_i)!=(0,0)
def wt(a,b):
    return sum(1 for i in range(n) if (a[i]%4,b[i]%4)!=(0,0))
# min nonzero weight in S (distance proxy lower bound for stabilizer elements)
mins=min(wt(a,b) for (a,b) in S if any(x for x in a) or any(x for x in b))
print("min weight nonzero stabilizer element:",mins)

# centralizer C: all (a|b) commuting with all 4 generators
def commutes_all(a,b):
    for i in range(4):
        if symp(np.array(a),np.array(b),A[i],B[i])!=0: return False
    return True
# enumerate full Pauli group 4^(2n)=4^10 ~ 1.05M -> feasible
C=[]
import numpy as np
count=0
minlog=99
# iterate
from itertools import product
for a in product(range(4),repeat=n):
    for b in product(range(4),repeat=n):
        if commutes_all(a,b):
            count+=1
            if (a,b) not in S:
                w=wt(a,b)
                if w<minlog: minlog=w
print("|C| =",count)
print("min weight of logical (C minus S) =",minlog)
# weight-1 logicals?
w1=0
for a in product(range(4),repeat=n):
    for b in product(range(4),repeat=n):
        if wt(a,b)==1 and commutes_all(a,b) and (a,b) not in S:
            w1+=1
print("number of weight-1 logicals:",w1)
