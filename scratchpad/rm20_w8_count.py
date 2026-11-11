"""Exact number of weight-8 words of RM_7(20,4) = RM(3,4)^perp on F_7^4 (8-point sets carrying a
full-support annihilator of the cubics), for the weight-8 mixed stage of [[1746,270,7]]_7.

Frame reduction as in rm20_w7_classify.py: up to AGL(4,7) every 8-set of affine span j contains the
frame 0, e_1..e_j plus 8-(j+1) further points of the coordinate j-flat (span 1 is impossible: a line
has 7 points).  After row-reducing the frame columns of the 35 x 2401 cubic evaluation matrix, the
extras' residuals decide dependence:
  span 4 (3 extras of 2396): weight 7 showed no zero residual and no proportional pair, so a set is
         dependent iff c' is a multiple of a' + beta b' (beta != 0): collisions between the
         6 C(2396,2) pair directions and the 2396 single directions;
  span 3 (4 extras of 339): weight 7 showed 60 dependent triples (line5+2) and no proportional pair;
         a 4-set is dependent iff it contains one of them or a' + beta b' ~ c' + delta d' for two
         disjoint pairs: collisions among the 6 C(339,2) pair directions;
  span 2 (5 extras of 46): batched rank of the 35 x 8 matrices (1.37e6 sets).
Each dependent set gets its null space and the number of full-support null vectors.

Orbit counting: an ordered affinely independent (j+1)-tuple f of a set S of span j fixes the affine
map onto the standard frame on span(S), so (S, f) <-> (R, h) with R a frame set and f = h(standard):
    #{S of span j with property X} = sum_{frame sets R with X} N_j / F(R),
N_j = #ordered affinely independent (j+1)-tuples of F_7^4, F(R) = #ordered affine frames inside R.
Checked on span 2: summing N_2/F(R) over all C(46,5) frame sets must give the number of 8-subsets
of the planes, 139 650 C(49,8).

Writes runs/weighted/rm20_w8_count.json.  usage: PYTHONPATH=. python scratchpad/rm20_w8_count.py
"""
import sys, os, json, time, itertools
from math import comb
from fractions import Fraction
import numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.pointsets import batched_rank_modp, frame_and_flat, inverse_table, direction_keys, _lexsort_rows
from qmsd.structured_ad import _null_space

p, m, W = 7, 4, 8
OUT = "runs/weighted/rm20_w8_count.json"
T0 = time.time()
def log(s): print(f"[{time.time()-T0:6.0f}s] {s}", flush=True)
state = json.load(open(OUT)) if os.path.exists(OUT) else {"p": p, "m": m, "weight": W, "spans": {}}
def save(): json.dump(state, open(OUT, "w"), indent=1)

P = np.asarray(points(m, p), dtype=np.int64)
R3 = np.asarray(rm_generator(3, m, p), dtype=np.int64) % p
assert R3.shape == (35, p ** m)
inv = inverse_table(p)
NPTS = p ** m

def N_tuples(j):
    """ordered affinely independent (j+1)-tuples of F_p^m"""
    v = NPTS
    for i in range(j):
        v *= NPTS - p ** i
    return v

def affine_rank(pts):
    pts = np.asarray(pts, dtype=np.int64)
    D = (pts[1:] - pts[0]) % p
    return int(batched_rank_modp(D[None], p)[0])

def frames_in(T, j):
    """number of ordered (j+1)-tuples of distinct points of T that are affinely independent"""
    pts = P[list(T)]
    tuples = np.array(list(itertools.permutations(range(len(T)), j + 1)), dtype=np.int64)
    D = (pts[tuples[:, 1:]] - pts[tuples[:, :1]]) % p          # (B, j, m)
    return int((batched_rank_modp(D, p) == j).sum())

def null_info(T):
    ns = np.asarray(_null_space(R3[:, list(T)], p)) % p
    k = ns.shape[0]
    if k == 0:
        return 0, 0
    assert k <= 6
    coefs = np.array(list(itertools.product(range(p), repeat=k)), dtype=np.int64)
    U = (coefs @ ns) % p
    return k, int((U != 0).all(axis=1).sum())

def frame_reduce(frame):
    M = R3.copy()
    for i, c in enumerate(frame):
        r = next(t for t in range(i, 35) if M[t, c] != 0)
        M[[i, r]] = M[[r, i]]
        M[i] = (M[i] * inv[M[i, c]]) % p
        for t in range(35):
            if t != i and M[t, c]:
                M[t] = (M[t] - M[t, c] * M[i]) % p
    return M

def collisions(keysA, keysB):
    """pairs (ia, ib) with keysA[ia] == keysB[ib] (rows)"""
    K = np.concatenate([keysA, keysB]); tag = np.concatenate([np.zeros(len(keysA), np.int64), np.ones(len(keysB), np.int64)])
    idx = np.concatenate([np.arange(len(keysA)), np.arange(len(keysB))])
    o = _lexsort_rows(K); K, tag, idx = K[o], tag[o], idx[o]
    same = np.ones(len(K), bool); same[1:] = (K[1:] == K[:-1]).all(axis=1); same[0] = False
    starts = np.nonzero(~same)[0]; ends = np.append(starts[1:], len(K))
    out = []
    for s, e in zip(starts, ends):
        if e - s < 2: continue
        a = idx[s:e][tag[s:e] == 0]; b = idx[s:e][tag[s:e] == 1]
        for x in a:
            for y in b:
                out.append((int(x), int(y)))
    return out

def pair_dirs(resid):
    n = len(resid)
    ii, jj = np.triu_indices(n, 1)
    V = np.concatenate([(resid[ii] + beta * resid[jj]) % p for beta in range(1, p)])
    I = np.tile(ii, p - 1); J = np.tile(jj, p - 1)
    nz = V.any(axis=1)
    return V[nz], I[nz], J[nz]

for j in (2, 3, 4):
    key = f"span{j}"
    if key in state["spans"] and state["spans"][key].get("done"):
        log(f"{key}: done, skipping"); continue
    t = time.time()
    frame, flat = frame_and_flat(m, j, p)
    others = np.array([x for x in flat if x not in set(frame)], dtype=np.int64)
    extra = W - (j + 1)
    M = frame_reduce(frame)
    resid = M[j + 1:, :][:, others].T.copy() % p
    assert resid.any(axis=1).all(), "zero residual"
    dep = set()
    rec = {"extras": extra, "others": len(others)}
    if j == 4:
        from qmsd.pointsets import _find_rows
        ks = direction_keys(resid, p); o = _lexsort_rows(ks); ks_sorted = ks[o]
        assert (ks_sorted[1:] != ks_sorted[:-1]).any(axis=1).all(), "proportional pair"
        n = len(resid); npairs = 0
        for x in range(n - 1):
            ys = np.arange(x + 1, n)
            for beta in range(1, p):
                V = (resid[x][None, :] + beta * resid[ys]) % p
                pos = _find_rows(ks_sorted, direction_keys(V, p))
                for yi in np.nonzero(pos >= 0)[0]:
                    c = int(o[pos[yi]]); y = int(ys[yi])
                    if c not in (x, y): dep.add(tuple(sorted((x, y, c))))
            npairs += len(ys)
            if x % 200 == 0: log(f"  {key}: x={x}/{n}, dependent so far {len(dep)}")
    elif j == 3:
        # dependent triples (all are line5+2 at weight 7) extended by any fourth point
        V, I, J = pair_dirs(resid)
        kp = direction_keys(V, p); ks = direction_keys(resid, p)
        triples = set()
        for a, c in collisions(kp, ks):
            x, y = int(I[a]), int(J[a])
            if c not in (x, y): triples.add(tuple(sorted((x, y, c))))
        rec["dependent_triples"] = len(triples)
        for tr in triples:
            for d in range(len(others)):
                if d not in tr: dep.add(tuple(sorted(tr + (d,))))
        # independent triples, dependent 4-set: a + beta b ~ c + delta d, disjoint pairs
        for a, b in collisions(kp, kp):
            if a >= b: continue
            s = {int(I[a]), int(J[a]), int(I[b]), int(J[b])}
            if len(s) == 4: dep.add(tuple(sorted(s)))
    else:
        B = 50_000; it = itertools.combinations(range(len(others)), extra)
        while True:
            ch = np.fromiter(itertools.chain.from_iterable(itertools.islice(it, B)), dtype=np.int64)
            if ch.size == 0: break
            ch = ch.reshape(-1, extra)
            rk = batched_rank_modp(resid[ch], p)
            for row in ch[rk < extra]: dep.add(tuple(int(v) for v in row))
    rec["dependent_sets"] = len(dep)
    log(f"{key}: {len(dep)} dependent sets [{time.time()-t:.0f}s]")
    # analyse
    Nj = N_tuples(j); total = Fraction(0); nfull = 0; hist = {}; examples = []
    for sset in sorted(dep):
        T = list(frame) + [int(others[i]) for i in sset]
        assert affine_rank(P[T]) == j
        k, full = null_info(T)
        assert k >= 1, "listed set is independent"
        if full == 0: continue
        F = frames_in(T, j)
        total += Fraction(Nj * full, F)
        nfull += 1
        lab = f"nullity{k}:full{full}:frames{F}"
        hist[lab] = hist.get(lab, 0) + 1
        if len(examples) < 10: examples.append([tuple(int(v) for v in P[x]) for x in T])
    assert total.denominator == 1, total
    rec.update({"full_support_sets": nfull, "types": hist, "examples": examples,
                "weight8_words": int(total.numerator), "N_j": Nj, "seconds": round(time.time() - t), "done": True})
    state["spans"][key] = rec; save()
    log(f"{key}: full-support sets {nfull}, types {hist}, weight-8 words {int(total)}")

# ---- check of the orbit-count formula on span 2: every 8-subset of a plane
if "orbit_check_span2" not in state:
    t = time.time()
    frame, flat = frame_and_flat(m, 2, p)
    others = [x for x in flat if x not in set(frame)]
    Ppl = P[flat]
    # collinear triples of the plane, as a set of sorted index triples (global indices)
    col = set()
    for a, b, c in itertools.combinations(flat, 3):
        if affine_rank(P[[a, b, c]]) < 2: col.add((a, b, c))
    acc = Fraction(0); Nj = N_tuples(2)
    for ext in itertools.combinations(others, 5):
        T = sorted(list(frame) + list(ext))
        nc = sum(1 for tr in itertools.combinations(T, 3) if tr in col)
        F = 8 * 7 * 6 - 6 * nc
        acc += Fraction(Nj, F)
    planes = 139650
    state["orbit_check_span2"] = {"sum": str(acc), "expected": planes * comb(49, 8), "ok": acc == planes * comb(49, 8),
                                  "seconds": round(time.time() - t)}
    save(); log(f"orbit check span 2: {state['orbit_check_span2']}")

tot = sum(v["weight8_words"] for v in state["spans"].values())
state["A8_RM20"] = tot; save()
log(f"A_8(RM_7(20,4)) = {tot}")
