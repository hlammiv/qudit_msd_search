"""p=11, m=4, r=13: hunt for a cubic g whose full zero set Z(g) is independent in degree 13
with |Z(g)| as close as possible to the cap 1319 = dim RM(13,4) - dim RM(10,4).

For each random cubic form g0 all eleven level sets {g0 = c} are cubic zero sets; the counts
are histogrammed and the ranks of RM(13,4) on the level sets with count in [lo, 1319] are
computed (galois).  Records every tested set to a JSONL file and the best one to a JSON.
Usage: python p11_cubic_hunt.py NCUBICS SEED LO OUT_PREFIX"""
import sys, time, json
import numpy as np
from qmsd.reedmuller import rm_generator, monomial_exponents
from hyp_common import rank_modp

p, m, r = 11, 4, 13
ncub = int(sys.argv[1]); seed = int(sys.argv[2]); lo = int(sys.argv[3]); pref = sys.argv[4]
rng = np.random.default_rng(seed)
G = np.asarray(rm_generator(r, m, p), dtype=np.int64)
G3 = np.asarray(rm_generator(3, m, p), dtype=np.int64)
exps3 = monomial_exponents(m, 3, p)
cap = G.shape[0] - rm_generator(10, m, p).shape[0]
print(f"cap = {cap}", flush=True)
best = None
counts_all = []
for t in range(ncub):
    c = rng.integers(0, p, size=len(exps3))
    top = [i for i, a in enumerate(exps3) if sum(a) == 3]
    if (c[top] == 0).all():
        c[top[0]] = 1
    vals = (c @ G3) % p
    cnt = np.bincount(vals, minlength=p)
    counts_all.extend(cnt.tolist())
    for lev in range(p):
        nZ = int(cnt[lev])
        if not (lo <= nZ <= cap):
            continue
        Z = np.nonzero(vals == lev)[0]
        t1 = time.time()
        h = rank_modp(G[:, Z], p)
        indep = (h == nZ)
        n = p ** m - nZ
        rec = dict(cubic=t, level=lev, nZ=nZ, h=h, independent=bool(indep), coeffs=c.tolist(),
                   gamma_d22=float(np.log(n / nZ) / np.log(22)), secs=round(time.time() - t1, 1))
        print(f"cubic {t} level {lev}: |Z|={nZ} h={h} {'INDEP' if indep else 'dep'} "
              f"gamma(d=22)={rec['gamma_d22']:.4f} [{rec['secs']}s]", flush=True)
        with open(pref + "_tested.jsonl", "a") as f:
            f.write(json.dumps(rec) + "\n")
        if indep and (best is None or nZ > best["nZ"]):
            best = dict(rec, puncture_columns_1indexed=(Z + 1).tolist(), p=p, m=m, r=r,
                        k=nZ, n=n, d_lower_bound=22,
                        proof="O = Z(g) for the cubic g (coeffs in monomial_exponents(4,3,11) order, "
                              "level set g = level); g*RM(10,4) lies in the check space, so every kernel "
                              "word e has g*e in RM(10,4)^perp = RM(29,4), whose minimum distance is 22.")
            with open(pref + "_best.json", "w") as f:
                json.dump(best, f)
counts_all = np.array(counts_all)
print("level-set counts: min", counts_all.min(), "max", counts_all.max(), "mean", counts_all.mean(),
      "fraction <= cap", (counts_all <= cap).mean(), "fraction in [lo,cap]", ((counts_all >= lo) & (counts_all <= cap)).mean())
if best:
    print("BEST:", best["k"], best["n"], best["gamma_d22"])
