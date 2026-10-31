"""Orbit-union output sets on the p = 7, m = 4, r = 7 sphere charge: the outputs are a union of
orbits of a finite group Gamma <= O(4, F_7) (matrices M with M^T M = I, so g(Mx) = g(x) and the
charge classes are preserved), of total size k, drawn at random from the cube-class points.
Cheap stages as in plane_sweep.py: plane-restricted distance d_plane (rigorous upper bound on
d), then the class-block weight <= 5 check for d_plane >= 7.  Survivors are candidates for the
full certificate (certify_candidate.py --file).  Appends JSON lines; resumable.

Why: a Gamma-symmetric code's low-weight words come in Gamma-orbits, so every existence scan
only needs leading columns among orbit representatives (|Gamma|-fold fewer), and the symmetric
sets are a different population from the random draws whose window closes at k = 270-275.

usage: orbit_sweep.py K SEED0 NSEEDS OUT.jsonl --group SPEC
  SPEC: '+'-joined generators from  shift (x1->x2->x3->x4), cyc3 (x1->x2->x3), swap12, swap34, neg (x -> -x), neg12,
        neg1 (x1 -> -x1), rot12 (rotation by (2,2) in the x1x2-plane), rot34, rot13.
        e.g. shift+neg (order 8), rot12+rot34 (order 64 -> use with care), shift (order 4)."""
import sys, json, random, itertools, time, os, numpy as np
from qmsd.reedmuller import points
from qmsd.weighted import weighted_support, quadric, class_points
from qmsd.structured_ad import _rank, _null_space
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import weight_d_exists_fast
k, seed0, nseeds, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
spec = sys.argv[sys.argv.index("--group") + 1] if "--group" in sys.argv else "shift+neg"
p, m, r = 7, 4, 7
P = np.asarray(points(m, p), dtype=np.int64)
index = {tuple(int(v) for v in x): i for i, x in enumerate(P)}
g = quadric(p, m, (1, 1, 1, 1)); GX, X, mu = weighted_support(p, m, g, r); cand = class_points(mu, p, 1)

# ---- the group
def gen(name):
    I = np.eye(m, dtype=np.int64)
    if name == "shift": return np.roll(I, 1, axis=0)
    if name == "swap12": M = I.copy(); M[[0, 1]] = M[[1, 0]]; return M
    if name == "swap34": M = I.copy(); M[[2, 3]] = M[[3, 2]]; return M
    if name == "neg": return (-I) % p
    if name == "neg1": M = I.copy(); M[0, 0] = p - 1; return M
    if name == "cyc3": return np.eye(m, dtype=np.int64)[[1, 2, 0, 3]]
    if name == "neg12": M = I.copy(); M[0, 0] = M[1, 1] = p - 1; return M
    if name.startswith("rot"):
        i, j = int(name[3]) - 1, int(name[4]) - 1; M = I.copy()
        M[i, i], M[i, j], M[j, i], M[j, j] = 2, p - 2, 2, 2        # (2,2): 2^2 + 2^2 = 8 = 1 mod 7
        return M
    raise ValueError(name)
gens = [gen(s) for s in spec.split("+")]
for M in gens: assert ((M.T @ M) % p == np.eye(m, dtype=np.int64)).all(), "generator is not orthogonal"
group = {tuple(np.eye(m, dtype=np.int64).flatten())}; frontier = [np.eye(m, dtype=np.int64)]
while frontier:
    new = []
    for A in frontier:
        for M in gens:
            B = (M @ A) % p; key = tuple(B.flatten())
            if key not in group: group.add(key); new.append(B)
    frontier = new
group = [np.array(t).reshape(m, m) for t in group]
# orbits on the class points
Xc = [int(X[l]) for l in cand]; cand_of_point = {Xc[i]: cand[i] for i in range(len(cand))}
seen = set(); orbits = []
for i, pt in enumerate(Xc):
    if pt in seen: continue
    orb = set()
    for M in group:
        q = index[tuple(int(v) for v in (M @ P[pt]) % p)]
        assert q in cand_of_point, "group element leaves the class (not orthogonal?)"
        orb.add(cand_of_point[q])
    seen.update(int(X[l]) for l in orb); orbits.append(sorted(orb))
sizes = sorted({len(o) for o in orbits})
print(f"group {spec}: order {len(group)}; {len(cand)} class points in {len(orbits)} orbits, sizes {sizes}", flush=True)

# ---- filters (as plane_sweep.py)
normals = np.array([v for v in itertools.product(range(p), repeat=m) if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1])
def rref2(a, b):
    M = np.array([a, b]) % p
    for c in range(m):
        piv = [i for i in range(2) if M[i, c]]
        if not piv: continue
        i = piv[0]
        if i: M[[0, i]] = M[[i, 0]]
        M[0] = (M[0] * pow(int(M[0, c]), -1, p)) % p; M[1] = (M[1] - M[1, c] * M[0]) % p; break
    for c in range(m):
        if M[1, c]:
            M[1] = (M[1] * pow(int(M[1, c]), -1, p)) % p; M[0] = (M[0] - M[0, c] * M[1]) % p; break
    return tuple(M.flatten())
subs = [np.array(s).reshape(2, m) for s in {rref2(normals[i], normals[j]) for i in range(len(normals)) for j in range(i + 1, len(normals))}]
PX = P[np.array(Xc)]
COEFS = {nu: np.array(list(itertools.product(range(p), repeat=nu)))[1:] for nu in range(1, 7)}
def min_null_weight(M, cap):
    ns = np.asarray(_null_space(M, p)) % p; nu = ns.shape[0]
    if nu == 0: return None
    if nu <= 6: return int(((COEFS[nu] @ ns) % p != 0).sum(axis=1).min())
    for dd in range(2, min(cap, 7)):
        if weight_d_exists_fast(M, p, dd): return dd
    return 7
def plane_distance(cols):
    G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p
    keep = [l for l in range(GX.shape[1]) if l not in set(cols)]; pos = {l: j for j, l in enumerate(keep)}
    sel = set(cols); surv = np.array([cand[i] not in sel for i in range(len(cand))])
    colidx = np.array([pos.get(cand[i], -1) for i in range(len(cand))])
    cls_cols = colidx[surv]; Gc = G0[:, cls_cols]; Gc = Gc[[t for t in range(Gc.shape[0]) if Gc[t].any()]]
    local = -np.ones(len(cand), dtype=int); local[np.nonzero(surv)[0]] = np.arange(surv.sum())
    best = 99
    for Sb in subs:
        v = (PX @ Sb.T) % p; key = v[:, 0] * p + v[:, 1]
        order = np.argsort(key, kind="stable"); ks = key[order]
        starts = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]]); ends = np.r_[starts[1:], len(ks)]
        for a, b in zip(starts, ends):
            idx = order[a:b]; cc = local[idx]; cc = cc[cc >= 0]
            if len(cc) < 2: continue
            M = Gc[:, cc]
            if _rank(M.T.copy(), p) == len(cc): continue
            w = min_null_weight(M, best)
            if w is not None and w < best:
                best = w
                if best < 7: return best, G0, Gc
    return best, G0, Gc

# ---- sampler: random union of orbits of total size exactly k, full rank
def draw(rng):
    for _ in range(200):
        order = orbits[:]; rng.shuffle(order); chosen = []; tot = 0
        for o in order:
            if tot + len(o) <= k: chosen.append(o); tot += len(o)
            if tot == k: break
        if tot != k: continue
        cols = sorted(c for o in chosen for c in o)
        if _rank(GX[:, cols].T.copy(), p) == k: return cols, sorted(len(o) for o in chosen)
    return None, None
done = set()
if os.path.exists(out):
    for line in open(out):
        if line.strip(): rec = json.loads(line); done.add((rec["k"], rec["seed"], rec.get("group")))
for seed in range(seed0, seed0 + nseeds):
    if (k, seed, spec) in done: continue
    rng = random.Random(1000 * k + seed)
    cols, osz = draw(rng)
    rec = {"k": k, "seed": seed, "group": spec, "group_order": len(group)}
    if cols is None:
        rec["d_plane"] = None; rec["note"] = "no full-rank orbit union of size k found"
    else:
        t0 = time.time(); dpl, G0, Gc = plane_distance(cols)
        rec.update({"n": G0.shape[1], "d_plane": dpl, "t_plane": round(time.time() - t0), "orbit_sizes": osz})
        if dpl >= 7:
            t1 = time.time(); db = None
            for dd in (2, 3, 4, 5):
                if weight_d_exists_fast(Gc, p, dd): db = dd; break
            rec["block_d_le5"] = db; rec["t_block"] = round(time.time() - t1); rec["cols_in_X"] = cols
    with open(out, "a") as f: f.write(json.dumps(rec) + "\n")
    print(json.dumps({kk: vv for kk, vv in rec.items() if kk != "cols_in_X"}), flush=True)
