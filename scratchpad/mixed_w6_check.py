"""Rigorous exclusion of weight-6 words that touch non-class columns, for a sphere-charge code.
(1) reproduce the k=260 sample-1 set from seed 7; (2) verify that every weight-5/6 word of
RM_7(20,4) (= RM(3,4)^perp) is supported on a line, by enumerating point sets up to affine
equivalence; (3) enumerate all line-supported weight-5/6 words with support in the non-class
points N, divide by w = (g-1)(g-6), and test whether any extends to a word of ker G0 by adding
at most one class column."""
import sys, json, random, itertools, time
import numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.weighted import weighted_support, quadric, class_points
from qmsd.structured_ad import _rank, _null_space
from qmsd.triorthogonal import shorten_matrix
p, m, r = 7, 4, 7
P = np.asarray(points(m, p), dtype=np.int64); index = {tuple(int(v) for v in x): i for i, x in enumerate(P)}
g = quadric(p, m, (1, 1, 1, 1)); GX, X, mu = weighted_support(p, m, g, r)
cand = class_points(mu, p, 1); cset = set(cand)
# ---- (1) reproduce the samples of wm4_bracket.py with seed 7, k list 275,265,260, 2 each
rng = random.Random(7); samples = []
for k in (275, 265, 260):
    for s in range(2):
        cols = sorted(rng.sample(cand, k))
        while _rank(GX[:, cols].T.copy(), p) < k:
            cols = sorted(rng.sample(cand, k))
        samples.append((k, s, cols))
ref = json.load(open("runs/weighted/p7_m4_sphere_bracket_2026-10-27.json"))
for (k, s, cols), rr in zip(samples, ref):
    assert (rr["k"], rr["sample"], rr["cols_in_X"]) == (k, s, cols), (k, s)
print("first 5 samples reproduce the lenore JSON; using k=260 sample 1", flush=True)
k, s, cols = samples[5]; assert (k, s) == (260, 1)
G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p
keep = [l for l in range(GX.shape[1]) if l not in set(cols)]            # G0 column -> index in X
n = G0.shape[1]; print("G0", G0.shape, flush=True)
json.dump({"charge": "sphere x1^2+x2^2+x3^2+x4^2", "p": p, "m": m, "r": r, "k": k, "n": n, "seed": "wm4_bracket.py seed 7, k=260 sample 1",
           "cols_in_X": cols, "X_global_1indexed": [int(x) + 1 for x in X], "mu": [int(u) for u in mu]},
          open("runs/weighted/p7_m4_sphere_1756_260_cols.json", "w"), indent=1)
# class / non-class columns of G0 (global point index)
gpt = [int(X[l]) for l in keep]                                            # global point index per G0 column
is_class = np.array([keep[j] in cset for j in range(n)])
C_cols = np.nonzero(is_class)[0]; N_cols = np.nonzero(~is_class)[0]
print("class cols", len(C_cols), "non-class cols", len(N_cols), flush=True)
# rows of G0 vanishing on the class columns = the kernel K; sanity: dim 35
K_dim = G0.shape[0] - _rank(G0[:, C_cols], p); print("dim of {f in G0: f|C = 0} =", K_dim, flush=True)
assert K_dim == 35
# ---- (2) every weight-5/6 word of RM(20,4) is supported on a line: frame enumeration
R3 = np.asarray(rm_generator(3, m, p), dtype=np.int64) % p                 # 35 x 2401 evaluation of RM(3,4)
assert R3.shape[0] == 35
def frame_pts(j):                                                          # 0, e1..ej as global indices
    return [index[(0, 0, 0, 0)]] + [index[tuple(1 if i == a else 0 for i in range(m))] for a in range(j)]
def flat_pts(j):                                                           # all points with coords j.. = 0
    return [i for i, x in enumerate(P) if all(int(v) == 0 for v in x[j:])]
t0 = time.time(); bad = 0; checked = 0; flagged = 0
for j in (2, 3, 4):
    fr = frame_pts(j); others = [i for i in flat_pts(j) if i not in set(fr)]
    for extra in (6 - (j + 1), 5 - (j + 1)):                               # |T| = 6 and 5
        if extra < 0: continue
        for comb in itertools.combinations(others, extra):
            T = fr + list(comb)
            checked += 1
            if _rank(R3[:, T].T.copy(), p) < len(T):
                ns = np.asarray(_null_space(R3[:, T], p)) % p          # does some null vector have FULL support?
                full = False
                for coef in itertools.product(range(p), repeat=ns.shape[0]):
                    if any(coef) and ((np.array(coef) @ ns) % p).all():
                        full = True; break
                flagged += 1
                if full: bad += 1
print(f"frame enumeration: {checked} point sets of span 2..4 checked, {flagged} dependent (proper-subset words), {bad} carry a full-support word [{time.time()-t0:.0f}s]", flush=True)
assert bad == 0, "lemma fails: non-collinear weight-5/6 words of RM(20,4) exist"
# ---- (3) line-supported weight-5/6 words of RM(20,4) with support in N
w = ((g - 1) * (g - 6)) % p                                                # charge weight per global point
Nset = set(gpt[j] for j in N_cols); col_of = {gpt[j]: j for j in range(n)}
winv = {x: pow(int(w[x]), -1, p) for x in Nset}
# directions of class columns for the weight-5 + one-class-column case
def norm(v):
    v = np.asarray(v) % p; nzv = np.nonzero(v)[0]
    if len(nzv) == 0: return None
    return tuple(((v * pow(int(v[nzv[0]]), -1, p)) % p).tolist())
cdir = {}
for j in C_cols:
    d = norm(G0[:, j]); assert d is not None and d not in cdir, "two class columns proportional"; cdir[d] = j
dirs = [v for v in itertools.product(range(p), repeat=m) if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1]
seen_lines = set(); n_lines = 0; cands = 0; hits = []
t0 = time.time()
for a in range(p ** m):
    for v in dirs:
        pts = tuple(sorted(index[tuple((int(P[a][i]) + t * v[i]) % p for i in range(m))] for t in range(p)))
        if pts in seen_lines: continue
        seen_lines.add(pts); n_lines += 1
        inN = [x for x in pts if x in Nset]
        if len(inN) < 5: continue
        ns = np.asarray(_null_space(R3[:, list(pts)], p)) % p               # annihilators of cubics on the 7 points
        assert ns.shape[0] == 3
        for coef in itertools.product(range(p), repeat=3):
            if not any(coef): continue
            u = (np.array(coef) @ ns) % p
            supp = [pts[i] for i in np.nonzero(u)[0]]
            if len(supp) not in (5, 6) or any(x not in Nset for x in supp): continue
            cands += 1
            eN = [(int(u[pts.index(x)]) * winv[x]) % p for x in supp]
            scols = [col_of[x] for x in supp]
            sv = (G0[:, scols] @ np.array(eN)) % p
            if len(supp) == 6:
                if not sv.any(): hits.append(("weight6_allN", scols, eN))
            else:
                d = norm((-sv) % p)
                if d is None: hits.append(("weight5_N_only?!", scols, eN))
                elif d in cdir: hits.append(("weight5_N+1class", scols, eN, int(cdir[d])))
print(f"lines {n_lines}, line words of weight 5/6 inside N: {cands}, extensions to ker G0: {len(hits)} [{time.time()-t0:.0f}s]", flush=True)
print("RESULT:", "no mixed weight-6 word" if not hits else hits[:5], flush=True)
json.dump({"k": k, "n": n, "frame_sets_checked": checked, "lines": n_lines, "line_word_candidates": cands, "mixed_hits": hits},
          open("runs/weighted/p7_m4_sphere_1756_260_mixed_w6_check.json", "w"), indent=1)
