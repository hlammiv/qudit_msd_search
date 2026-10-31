"""Which symmetry groups allow a full-rank orbit-union output set of size k on the sphere charge?
(Orbit unions respect the isotypic decomposition of RM(7,4)|_X under Gamma, so a component that is
too small caps the rank: e.g. -I in Gamma forces k <= 2 dim(even part).)  Reports, per group, the
rank deficiency of random unions at the given k's."""
import sys, random, numpy as np, importlib.util
from qmsd.reedmuller import points
from qmsd.weighted import weighted_support, quadric, class_points
from qmsd.structured_ad import _rank
p, m, r = 7, 4, 7
P = np.asarray(points(m, p), dtype=np.int64); index = {tuple(int(v) for v in x): i for i, x in enumerate(P)}
GX, X, mu = weighted_support(p, m, quadric(p, m, (1, 1, 1, 1)), r); cand = class_points(mu, p, 1)
def gen(name):
    I = np.eye(m, dtype=np.int64)
    if name == "shift": return np.roll(I, 1, axis=0)
    if name == "swap12": M = I.copy(); M[[0, 1]] = M[[1, 0]]; return M
    if name == "swap34": M = I.copy(); M[[2, 3]] = M[[3, 2]]; return M
    if name == "neg": return (-I) % p
    if name == "neg1": M = I.copy(); M[0, 0] = p - 1; return M
    if name == "neg12": M = I.copy(); M[0, 0] = M[1, 1] = p - 1; return M
    if name.startswith("rot"):
        i, j = int(name[3]) - 1, int(name[4]) - 1; M = I.copy(); M[i, i], M[i, j], M[j, i], M[j, j] = 2, p - 2, 2, 2; return M
    if name == "cyc3": M = np.eye(m, dtype=np.int64)[[1, 2, 0, 3]]; return M
    raise ValueError(name)
def closure(gens):
    group = {tuple(np.eye(m, dtype=np.int64).flatten())}; frontier = [np.eye(m, dtype=np.int64)]
    while frontier:
        new = []
        for A in frontier:
            for M in gens:
                B = (M @ A) % p; key = tuple(B.flatten())
                if key not in group: group.add(key); new.append(B)
        frontier = new
    return [np.array(t).reshape(m, m) for t in group]
Xc = [int(X[l]) for l in cand]; cand_of_point = {Xc[i]: cand[i] for i in range(len(cand))}
def orbits_of(group):
    seen = set(); orbs = []
    for pt in Xc:
        if pt in seen: continue
        orb = {cand_of_point[index[tuple(int(v) for v in (M @ P[pt]) % p)]] for M in group}
        seen.update(int(X[l]) for l in orb); orbs.append(sorted(orb))
    return orbs
ks = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [272, 280]
for spec in ["shift", "cyc3", "swap12", "swap12+swap34", "shift+swap12", "neg1", "neg12", "neg1+swap12", "rot12", "rot12+rot34", "rot12+shift", "shift+neg"]:
    group = closure([gen(s) for s in spec.split("+")]); orbs = orbits_of(group)
    sizes = sorted({len(o) for o in orbs}); rng = random.Random(1)
    line = f"{spec:14s} |G|={len(group):3d} orbits={len(orbs):3d} sizes={sizes}"
    for k in ks:
        defs = []
        for t in range(6):
            order = orbs[:]; rng.shuffle(order); chosen = []; tot = 0
            for o in order:
                if tot + len(o) <= k: chosen.append(o); tot += len(o)
                if tot == k: break
            if tot != k: continue
            cols = sorted(c for o in chosen for c in o); defs.append(k - _rank(GX[:, cols].T.copy(), p))
        line += f" | k={k}: deficiency {defs}"
    print(line, flush=True)
