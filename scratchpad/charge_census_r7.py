"""Degree-2 charges at r = 7 (p = 7, m = 4): for each class (unit 1, 2, 4) the class-point count
N_c, the rank rho of RM(7,4)|X on the class points (k ceiling; block checks = rho - k), and the
accidental-word closure k* = largest k with E[weight-6 block words] = C(N_c-k,6) 6^5 / 7^(rho-k)
< 0.1, with gamma at k* for d = 7 (n = |X| - k*)."""
import itertools, numpy as np
from math import comb, log
from qmsd.weighted import weighted_support, quadric, class_points
from qmsd.structured_ad import _rank
p, m, r = 7, 4, 7
def kstar(Nc, rho, X):
    best = None
    for k in range(100, rho):
        if comb(Nc - k, 6) * 6 ** 5 / 7 ** (rho - k) < 0.1: best = k
    return best
rows = []
forms = []
for diag in [(1, 1, 1, 1), (1, 1, 1, 3), (1, 1, 3, 3), (1, 3, 3, 3), (3, 3, 3, 3), (1, 1, 1, 0), (1, 1, 3, 0), (1, 3, 3, 0), (1, 1, 0, 0), (1, 3, 0, 0)]:
    for const in range(7):
        forms.append((f"{diag}+{const}", quadric(p, m, diag, const=const)))
    if 0 in diag:
        forms.append((f"{diag}+x{diag.index(0)+1}", quadric(p, m, diag, lin=tuple(1 if d == 0 else 0 for d in diag))))
print(f"{'charge':18s} |X|  rank | unit  N_c  rho  vanish  k*   n/k   gamma(d=7)")
for name, g in forms:
    GX, X, mu = weighted_support(p, m, g, r); rk = _rank(GX, p)
    for u in (1, 2, 4):
        cls = class_points(mu, p, u)
        if len(cls) < 300: continue
        rho = _rank(GX[:, cls], p); ks = kstar(len(cls), rho, len(X))
        if ks is None: continue
        gam = log((len(X) - ks) / ks) / log(7)
        print(f"{name:18s} {len(X):4d} {rk:4d} | {u}   {len(cls):4d} {rho:4d} {rk-rho:5d}  {ks:4d} {(len(X)-ks)/ks:5.2f}  {gam:.3f}", flush=True)
