"""Mixed weight-7 stage for a sphere-charge candidate (p = 7, m = 4, r = 7): exclude every weight-7
word of ker G0 that touches a non-class column.

For e in ker G0 the 35 kernel rows (g-1)(g-6) RM(3,4) force u = w * e_N (w = (g-1)(g-6), e_N the
non-class part) into RM(3,4)^perp = RM_7(20,4), whose words of weight <= 7 are supported on
lines (weight 5, 6: mixed_w6_check.py; weight 7: rm20_w7_classify.py -- the only 7-point sets
carrying a full-support annihilator of the cubics are full lines).  On a line the dual words
are the 57 projective words of RM(2,1): 21 of weight 5 (two roots), 14 of weight 6 (one root:
linear, or a square), 22 of weight 7 (constants and the 21 irreducible monic quadratics).

So a weight-7 word with e_N != 0 is one of
    (a) a weight-7 line word with all 7 points in N, no class column            (s = 0),
    (b) a weight-6 line word with 6 points in N + one class column               (s ~ c_i),
    (c) a weight-5 line word with 5 points in N + two class columns              (s in <c_i, c_j>),
where s is the residual syndrome G0 e_N projected to the 21 block coordinates (the coordinates
in which the class columns live after the kernel rows are split off).  (a)-(c) are one hash
lookup each: the directions of all single class columns and of all 6 (p-1) combinations
c_i + beta c_j of each pair are tabulated once (402 + 6 C(402,2) = 484 008 keys in F_7^21), and
every candidate's normalised s is looked up by binary search.  Lower-weight hits (s = 0 for a
weight-5/6 word, s ~ c_i for a weight-5 word) are reported too; they would contradict the
mixed6 stage.

usage: PYTHONPATH=. python scratchpad/mixed_w7_check.py K SEED [OUT.json] [--cols FILE.json]
       (the draw is read from runs/weighted/plane_sweep/plane_sweep_k*.jsonl as in certify_candidate.py,
        or from --cols FILE.json with a cols_in_X list, e.g. runs/weighted/p7_m4_sphere_1756_260_cols.json)
"""
import sys, os, json, time, itertools, glob
import numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.weighted import weighted_support, quadric, class_points
from qmsd.structured_ad import _rank, _null_space, _rref
from qmsd.triorthogonal import shorten_matrix
from qmsd.pointsets import affine_lines, inverse_table, pair_direction_table, lookup_directions

k, seed = int(sys.argv[1]), int(sys.argv[2])
args = [a for a in sys.argv[3:] if not a.startswith("--")]
out = args[0] if args else f"runs/certify/p7_k{k}_seed{seed}_mixed7.json"
cols_file = sys.argv[sys.argv.index("--cols") + 1] if "--cols" in sys.argv else None     # JSON with cols_in_X (draws outside the sweep)
p, m, r = 7, 4, 7
T0 = time.time()
def log(msg): print(f"[{time.time()-T0:6.0f}s] {msg}", flush=True)
inv = inverse_table(p)

# ---- the draw (same loading as certify_candidate.py)
rec = None
if cols_file:
    rec = json.load(open(cols_file)); assert rec["k"] == k
else:
    for f in glob.glob("runs/weighted/plane_sweep/plane_sweep_k*.jsonl"):
        for l in open(f):
            if not l.strip(): continue
            rr = json.loads(l)
            if rr["k"] == k and rr["seed"] == seed and "cols_in_X" in rr: rec = rr
    assert rec is not None, "survivor not found in the sweep files"
cols = rec["cols_in_X"]
cert = f"runs/certify/p7_k{k}_seed{seed}.json"
if os.path.exists(cert) and not cols_file:
    assert json.load(open(cert))["cols_in_X"] == cols, "sweep record and certificate disagree on the draw"
P = np.asarray(points(m, p), dtype=np.int64)
g = quadric(p, m, (1, 1, 1, 1)); GX, X, mu = weighted_support(p, m, g, r)
cand = class_points(mu, p, 1); cset = set(cand)
G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p
keep = [l for l in range(GX.shape[1]) if l not in set(cols)]
n = G0.shape[1]
is_class = np.array([keep[j] in cset for j in range(n)])
C_cols = np.nonzero(is_class)[0]; N_cols = np.nonzero(~is_class)[0]
gpt = np.array([int(X[l]) for l in keep])                                 # G0 column -> global point
log(f"k={k} seed={seed}: G0 {G0.shape}, class cols {len(C_cols)}, non-class cols {len(N_cols)}")

# ---- split the 56 rows into 21 block coordinates + 35 kernel rows: rref of [Gc | I]
Gc = G0[:, C_cols]
R, piv = _rref(np.concatenate([Gc, np.eye(G0.shape[0], dtype=np.int64)], axis=1), p)
rb = sum(1 for c in piv if c < Gc.shape[1])                              # rank of the class block
Mrow = R[:, Gc.shape[1]:]                                                 # invertible: M @ Gc = R[:, :402]
G0p = (Mrow @ G0) % p
assert (G0p[rb:, :][:, C_cols] == 0).all() and _rank(G0p[rb:], p) == G0.shape[0] - rb
K_dim = G0.shape[0] - rb
log(f"class block rank {rb}, kernel rows {K_dim} (expect 35)")
assert K_dim == 35
Cb = G0p[:rb, :][:, C_cols]                                               # 21 x 402 block coordinates of the class columns
assert (Cb != 0).any(axis=0).all()

# ---- direction tables: zero, single class columns, pairs c_i + beta c_j (qmsd.pointsets)
table = pair_direction_table(Cb, p)
skeys = table["keys"]
dups = int((skeys[1:] == skeys[:-1]).all(axis=1).sum())
log(f"direction table: {len(C_cols)} singles + {(p - 1) * len(C_cols) * (len(C_cols) - 1) // 2} pair directions (+ zero) = {len(skeys)} keys, duplicates {dups}")

# ---- line words: the 57 projective words of RM(2,1) on t = 0..6
R3 = np.asarray(rm_generator(3, m, p), dtype=np.int64) % p
ns = np.asarray(_null_space(R3[:, :p], p)) % p                            # points (t,0,0,0) have index t
assert ns.shape[0] == 3
Q = []; seen = set()
for coef in itertools.product(range(p), repeat=3):
    if not any(coef): continue
    u = (np.array(coef) @ ns) % p; f = np.nonzero(u)[0][0]; u = (u * inv[u[f]]) % p
    if tuple(u) in seen: continue
    seen.add(tuple(u)); Q.append(u)
Q = np.array(Q); Qw = (Q != 0).sum(axis=1)
assert len(Q) == 57 and dict(zip(*np.unique(Qw, return_counts=True))) == {5: 21, 6: 14, 7: 22}

# ---- non-class points, their columns and charge weights
w = ((g - 1) * (g - 6)) % p
inN = np.zeros(p ** m, dtype=bool); inN[gpt[N_cols]] = True
col_of = -np.ones(p ** m, dtype=np.int64); col_of[gpt] = np.arange(n)
winv = np.zeros(p ** m, dtype=np.int64); winv[inN] = inv[w[inN]]
LINES = affine_lines(p, m)                                                # (137200, 7), parameter order
nN = inN[LINES].sum(axis=1)
log(f"lines: {len(LINES)}; with 7/6/5 points in N: {int((nN == 7).sum())}/{int((nN == 6).sum())}/{int((nN == 5).sum())}")

# ---- candidates and the lookup
BATCH = 20_000
counts = {5: 0, 6: 0, 7: 0}; hits = []
tested_hits = 0; planted = []
for qi in range(len(Q)):
    q = Q[qi]; S = np.nonzero(q)[0]; wq = int(Qw[qi])
    lines = np.nonzero(inN[LINES[:, S]].all(axis=1))[0]
    counts[wq] += len(lines)
    for lo in range(0, len(lines), BATCH):
        L = LINES[lines[lo:lo + BATCH]][:, S]                              # (B, wq) global points
        e = (q[S][None, :] * winv[L]) % p                                  # (B, wq) coefficients e_N = u / w
        cols_b = col_of[L]
        assert (cols_b >= 0).all()
        syn = np.einsum("bwr,bw->br", G0p.T[cols_b], e) % p                # (B, 56)
        assert not syn[:, rb:].any(), "kernel rows do not vanish: not an RM(20,4) word"
        sb = syn[:, :rb]
        if lo == 0:                                                        # serial recheck of the batched syndromes
            for bi in np.random.default_rng(qi).integers(0, len(L), 3):
                assert ((G0p[:, cols_b[bi]] @ e[bi]) % p == syn[bi]).all()
            planted.append(sb[0].copy())
        pos = lookup_directions(table, sb, p)                               # zero rows -> the zero key
        for bi in np.nonzero(pos >= 0)[0]:
            oi = pos[bi]
            # verify explicitly: solve e_N + class part
            vec = np.zeros(n, dtype=np.int64); vec[cols_b[bi]] = e[bi]
            extra = [int(c) for c in (table["i"][oi], table["j"][oi]) if c >= 0]   # indices into C_cols
            if extra:
                A = np.concatenate([sb[bi][:, None], Cb[:, extra]], axis=1)
                sol = np.asarray(_null_space(A, p)) % p
                sol = sol[np.nonzero(sol[:, 0])[0][0]]; sol = (sol * inv[sol[0]]) % p
                for c, a in zip(extra, sol[1:]): vec[C_cols[c]] = a
                extra = [int(C_cols[c]) for c in extra]
            assert ((G0 @ vec) % p == 0).all(), "lookup hit does not verify"
            wt = int((vec != 0).sum()); tested_hits += 1
            hits.append({"weight": wt, "line": int(lines[lo + bi]), "line_word_weight": wq, "class_cols": extra,
                         "G0_cols": [int(c) for c in np.nonzero(vec)[0]], "coeffs": [int(vec[c]) for c in np.nonzero(vec)[0]]})
    if qi % 10 == 9: log(f"  {qi + 1}/57 line words done; candidates so far {counts}, hits {len(hits)}")
total = sum(counts.values())
log(f"candidates: weight-5 words in N {counts[5]}, weight-6 {counts[6]}, weight-7 {counts[7]} (total {total}); hits {len(hits)}")
# ---- positive control: plant the syndromes of real candidates as fake class columns and look them up
rng = np.random.default_rng(1)
planted = [s for s in planted if s.any()][:4]
assert len(planted) >= 3
ca = rng.integers(0, p, rb)
Cfake = np.concatenate([Cb, (2 * planted[0])[:, None], ca[:, None], ((planted[1] - 3 * ca) % p)[:, None]], axis=1)
tab2 = pair_direction_table(Cfake, p)
pos2 = lookup_directions(tab2, np.stack(planted[:3]), p)
nC = len(C_cols)
assert tab2["kind"][pos2[0]] == 1 and tab2["i"][pos2[0]] == nC, "planted single column not found"
assert tab2["kind"][pos2[1]] == 2 and {int(tab2["i"][pos2[1]]), int(tab2["j"][pos2[1]])} == {nC + 1, nC + 2}, "planted pair not found"
assert lookup_directions(table, np.stack(planted[:3]), p).max() < 0
log("positive control: planted single / pair class columns are found by the lookup")
verdict = "no mixed word of weight <= 7" if not hits else f"{len(hits)} mixed words, min weight {min(h['weight'] for h in hits)}"
log("RESULT: " + verdict)
json.dump({"k": k, "seed": seed, "n": n, "class_cols": int(len(C_cols)), "nonclass_cols": int(len(N_cols)),
           "block_rank": int(rb), "kernel_dim": int(K_dim), "direction_keys": int(len(skeys)), "duplicate_keys": dups,
           "lines_with_7_6_5_in_N": [int((nN == 7).sum()), int((nN == 6).sum()), int((nN == 5).sum())],
           "candidates_by_line_word_weight": {str(a): b for a, b in counts.items()}, "hits": hits,
           "verdict": verdict, "seconds": round(time.time() - T0)}, open(out, "w"), indent=1)
log("saved " + out)
