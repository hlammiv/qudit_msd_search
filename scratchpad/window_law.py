"""The accidental-word law for weighted sphere-charge supports: for (p, m, r) compute |X|, the rank,
the output class (cube class when 3 | p-1, else everything), the class-restricted rank, and the
predicted window: E[projective weight-w words on the class block] = C(N_class - k, w)(p-1)^(w-1)/p^(rank_class - k).
usage: window_law.py p m r"""
import sys, time, numpy as np
from math import comb, exp, log
from qmsd.weighted import weighted_support, quadric, class_points
from qmsd.structured_ad import _rank
p, m, r = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); t0 = time.time()
GX, X, mu = weighted_support(p, m, quadric(p, m, tuple([1] * m)), r)
cls = class_points(mu, p, 1) if (p - 1) % 3 == 0 else list(range(len(X)))
rk = _rank(GX.copy(), p); rk_cls = _rank(np.ascontiguousarray(GX[:, cls]).copy(), p) if len(cls) < len(X) else rk
print(f"p={p} m={m} r={r}: |X|={len(X)} rank={rk} class points={len(cls)} class rank={rk_cls} kernel rows={rk-rk_cls} ({time.time()-t0:.0f}s)", flush=True)
print("   k      n  classcols checks   E[w5]    E[w6]   P(no w6)   E[w7]   E[w8]  g(d=6)  g(d=7)  g(d=8)")
def E(N, w, chk): return comb(N, w) * (p - 1) ** (w - 1) / p ** chk
# scan k from where E[w6] ~ 1e-4 up to where E[w5] ~ 1
ks = [k for k in range(1, rk_cls) if 1e-5 < E(len(cls) - k, 6, rk_cls - k) < 1e3 or 1e-3 < E(len(cls) - k, 5, rk_cls - k) < 10]
for k in ks[:: max(1, len(ks) // 14)]:
    chk = rk_cls - k; N = len(cls) - k; n = len(X) - k
    print(f"{k:5d} {n:7d} {N:9d} {chk:6d}  {E(N,5,chk):8.1e} {E(N,6,chk):8.1e}  {exp(-E(N,6,chk)):6.2f}  {E(N,7,chk):8.1e} {E(N,8,chk):8.1e}  {log(n/k)/log(6):.3f}  {log(n/k)/log(7):.3f}  {log(n/k)/log(8):.3f}", flush=True)
