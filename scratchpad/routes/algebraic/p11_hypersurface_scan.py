"""p=11, m=4, r=13: sample polynomials g of degree delta, compute Z(g) and the rank h of
RM(13,4) on Z(g).  If h == |Z|, O = Z(g) is a valid output set with provable
d >= d(RM(26+delta,4)) (delta=3: 22, delta=4: 11, delta=5: 10, ...).
Usage: python p11_hypersurface_scan.py DELTA TRIALS SEED [outfile.jsonl]"""
import sys, time, json
import numpy as np
from qmsd.reedmuller import rm_generator, d_rm
from hyp_common import random_poly_values, rank_modp, bound

p, m, r = 11, 4, 13
delta = int(sys.argv[1]); trials = int(sys.argv[2]); seed = int(sys.argv[3])
out = sys.argv[4] if len(sys.argv) > 4 else None
rng = np.random.default_rng(seed)
t0 = time.time()
G = np.asarray(rm_generator(r, m, p), dtype=np.int64)
K = G.shape[0]
cap = K - rm_generator(r - delta, m, p).shape[0]
b = bound(delta, r, m, p)
print(f"dim RM(13,4)={K}, deg g={delta}: cap on k = {cap}, theorem bound d >= {b}; setup {time.time()-t0:.1f}s", flush=True)
for t in range(trials):
    g, c = random_poly_values(delta, m, p, rng, homogeneous_top=True)
    Z = np.nonzero(g == 0)[0]
    t1 = time.time()
    h = rank_modp(G[:, Z], p)
    n = p**m - len(Z)
    rec = dict(delta=delta, seed=seed, trial=t, nZ=int(len(Z)), h=int(h), cap=cap, bound=b,
               independent=bool(h == len(Z)), coeffs=c.tolist(),
               gamma_if_full=(float(np.log(n / len(Z)) / np.log(b)) if h == len(Z) else None))
    print(f"trial {t}: |Z|={len(Z)} h_Z(13)={h} {'INDEPENDENT k=%d gamma<=%.3f' % (len(Z), rec['gamma_if_full']) if rec['independent'] else 'dependent (excess %d)' % (len(Z)-h)}  [{time.time()-t1:.1f}s]", flush=True)
    if out:
        with open(out, "a") as f: f.write(json.dumps(rec) + "\n")
