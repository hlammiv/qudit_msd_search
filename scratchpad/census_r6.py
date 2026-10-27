"""Degree-5 charges at r = 6, p = 7, m = 4 (budget 23 - 18 = 5): support size, rank of RM(6,4) on
the support, cube-class points and the rank on them (the k ceiling), and the dimension of the
stabilizers vanishing on the whole class (the class-block weakness)."""
import itertools, numpy as np
from qmsd.reedmuller import points
from qmsd.weighted import weighted_support, quadric, product_of_linear, charge_degree_budget
from qmsd.structured_ad import _rank
from qmsd.triorthogonal import is_cube_mod_p
p, m, r = 7, 4, 6
print("budget deg g <=", charge_degree_budget(p, m, r), flush=True)
P = np.asarray(points(m, p), dtype=np.int64)
e = lambda i: tuple(1 if j == i else 0 for j in range(m))
forms = {
  "sphere x.x (deg 2)":               quadric(p, m, (1, 1, 1, 1)),
  "x1(x1-1)...(x1-4)":                product_of_linear(p, m, [(e(0), c) for c in range(5)]),
  "x1 x2 x3 x4 (x1+x2+x3+x4-1)":      product_of_linear(p, m, [(e(0), 0), (e(1), 0), (e(2), 0), (e(3), 0), ((1, 1, 1, 1), 1)]),
  "x1 x2 x3 x4 (x1+2x2+3x3+4x4-1)":   product_of_linear(p, m, [(e(0), 0), (e(1), 0), (e(2), 0), (e(3), 0), ((1, 2, 3, 4), 1)]),
  "x1(x1-1) x2(x2-1) x3":             product_of_linear(p, m, [(e(0), 0), (e(0), 1), (e(1), 0), (e(1), 1), (e(2), 0)]),
  "x1(x1-1)(x1-2) x2 x3":             product_of_linear(p, m, [(e(0), 0), (e(0), 1), (e(0), 2), (e(1), 0), (e(2), 0)]),
  "5 general hyperplanes":            product_of_linear(p, m, [((1, 0, 0, 0), 0), ((0, 1, 0, 0), 0), ((0, 0, 1, 0), 0), ((1, 1, 1, 1), 2), ((1, 2, 4, 1), 3)]),
  "sphere * x1 x2 x3":                (quadric(p, m, (1, 1, 1, 1)) * product_of_linear(p, m, [(e(0), 0), (e(1), 0), (e(2), 0)])) % p,
  "(sphere+1) * x1(x1-1)(x1-2)":      ((quadric(p, m, (1, 1, 1, 1), const=1)) * product_of_linear(p, m, [(e(0), 0), (e(0), 1), (e(0), 2)])) % p,
  "sphere * (sphere+1) * x1":         (quadric(p, m, (1, 1, 1, 1)) * quadric(p, m, (1, 1, 1, 1), const=1) * product_of_linear(p, m, [(e(0), 0)])) % p,
  "sum x_i^5 + 1":                    (np.sum(P ** 5, axis=1) + 1) % p,
  "x1^5 + x2^3 + x3^2 + x4":          (P[:, 0] ** 5 + P[:, 1] ** 3 + P[:, 2] ** 2 + P[:, 3]) % p,
  "x1^4 x2 + x3^2 x4 + 1":            (P[:, 0] ** 4 * P[:, 1] + P[:, 2] ** 2 * P[:, 3] + 1) % p,
}
print(f"{'g':34s} |X|  rank | cube pts  rank_on_cube  vanish_on_cube | 2cube  4cube")
for name, g in forms.items():
    GX, X, mu = weighted_support(p, m, g, r)
    rk = _rank(GX, p)
    cls = {u: [l for l, v in enumerate(mu) if is_cube_mod_p(int(v) * pow(u, -1, p) % p, p)] for u in (1, 2, 4)}
    rc = _rank(GX[:, cls[1]], p)
    print(f"{name:34s} {len(X):4d} {rk:4d} | {len(cls[1]):5d} {rc:9d} {rk - rc:12d} | {len(cls[2]):5d} {len(cls[4]):5d}", flush=True)
