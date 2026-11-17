"""Weight-8 census of a paired-cap code [[686 - 2t, 2t, 8]]_7 (probe_paired_cap.py, d' = 4, d'' >= 8).

By d = min(2d', d'') every weight-8 word of ker(G0) is one of
  (i)  (a', b') with a' on x1 = 5 and b' on x1 = 6 both weight-4 words of ker(RM(5,3) shortened at Y) ("d'-words"), and
       G0 (a', 0) + G0 (0, b') = 0; this includes the paired (a' | -a');
  (ii) a weight-8 word on one hyperplane only, i.e. a weight-8 word of ker(RM(6,3) shortened at Y) ("d''-words").
(i) is exact here: all d'-words by a 2+2 meet-in-the-middle, then syndrome matching in G0.
(ii) structural family only: off-Y parts of two lines carrying 6 points of Y between them (each such 8-set is tested for a
full-support dependency).  Completeness of (ii) is NOT claimed; see isd_paired_w8.py.
Every word is classified against the canonical logical rows H1: per-output failure counts give L_8.
usage: PYTHONPATH=. python scratchpad/census_paired_w8.py SWEEP.json T SAMPLE OUT.json [EXACT_II.json from w8_mitm44v2.py]"""
import sys, json, time, itertools, numpy as np
from collections import defaultdict
from qmsd.reedmuller import points, rm_generator
from qmsd.weighted import weighted_support, product_of_linear, weighted_code
from qmsd.triorthogonal import shorten_matrix
from qmsd.structured_ad import _rref, _null_space

src, t, s, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]; p = 7
T0 = time.time()
dr = [d for d in json.load(open(src))["draws"] if d["t"] == t and d["sample"] == s][0]; Y = dr["Y"]; Yset = set(Y)
P3 = np.asarray(points(3, p), dtype=np.int64); idx3 = {tuple(int(v) for v in x): i for i, x in enumerate(P3)}
offY = [x for x in range(p ** 3) if x not in Yset]; posY = {x: j for j, x in enumerate(offY)}

def basis(M):
    R, piv = _rref(np.asarray(M, dtype=np.int64) % p, p); return R[:len(piv)]

R5 = np.asarray(rm_generator(5, 3, p), dtype=np.int64) % p; R6 = np.asarray(rm_generator(6, 3, p), dtype=np.int64) % p
H5 = basis(shorten_matrix(R5, Y, p)); H6 = basis(shorten_matrix(R6, Y, p))
print(f"t={t}: off-Y columns {len(offY)}, rank H5 {H5.shape[0]}, rank H6 {H6.shape[0]}", flush=True)

def code(v):                                       # syndrome vector -> int (base p)
    return int(sum(int(c) * p ** i for i, c in enumerate(v)))

# ---- d'-words: all weight-4 dependencies of H5 (projective, first coefficient 1), by 2+2 MITM
N = len(offY); r5 = H5.shape[0]; W = p ** np.arange(r5, dtype=np.int64)
pairs = np.array(list(itertools.combinations(range(N), 2)), dtype=np.int64)
Hc = H5.T % p                                     # N x r5
left = {}                                          # code -> list of (i, j, cj) with coeffs (1, cj)
for c in range(1, p):
    syn = (Hc[pairs[:, 0]] + c * Hc[pairs[:, 1]]) % p
    for (i, j), key in zip(pairs, syn @ W):
        left.setdefault(int(key), []).append((int(i), int(j), c))
words = set()
for a in range(1, p):
    for b in range(1, p):
        syn = (-(a * Hc[pairs[:, 0]] + b * Hc[pairs[:, 1]])) % p
        keys = syn @ W
        for (k, l), key in zip(pairs, keys):
            for (i, j, cj) in left.get(int(key), ()):
                if len({i, j, int(k), int(l)}) < 4: continue
                sup = {i: 1, j: cj, int(k): a, int(l): b}
                cols = sorted(sup); first = sup[cols[0]]; inv = pow(first, p - 2, p)
                words.add(tuple((cc, sup[cc] * inv % p) for cc in cols))
words = sorted(words)
assert all(not (H5[:, [c for c, _ in w]] @ np.array([v for _, v in w]) % p).any() for w in words)
print(f"d'-words (weight 4, projective): {len(words)} [{time.time()-T0:.0f}s]", flush=True)

# ---- full code
g = product_of_linear(p, 4, [((1, 0, 0, 0), c) for c in range(5)]); GX, X, mu = weighted_support(p, 4, g, 6)
P4 = np.asarray(points(4, p), dtype=np.int64)[X]; idxX = {tuple(int(v) for v in x): i for i, x in enumerate(P4)}
S = sorted(idxX[(c,) + tuple(int(v) for v in P3[y])] for y in Y for c in (5, 6))
G0, H1, coeffs = weighted_code(GX, mu, S, p)
keep = [l for l in range(len(X)) if l not in set(S)]; kpos = {l: j for j, l in enumerate(keep)}
G0b = basis(G0); n = len(keep); kk = H1.shape[0]
print(f"full code: n={n} k={kk} rank G0={G0b.shape[0]}", flush=True)

def full_col(hyper, c3):                           # off-Y column c3 on hyperplane x1 = hyper -> column of the full code
    return kpos[idxX[(hyper,) + tuple(int(v) for v in P3[offY[c3]])]]

def vec(entries):
    e = np.zeros(n, dtype=np.int64)
    for col, val in entries: e[col] = (e[col] + val) % p
    return e

def proj_key(v):
    nz = np.nonzero(v)[0]
    if len(nz) == 0: return (0,)
    inv = pow(int(v[nz[0]]), p - 2, p); return tuple(int(x) for x in (v * inv) % p)

# (i): match G0-syndromes of the lifts, projectively
syn5 = [(G0b @ vec([(full_col(5, c), a) for c, a in w])) % p for w in words]
syn6 = [(G0b @ vec([(full_col(6, c), a) for c, a in w])) % p for w in words]
by6 = defaultdict(list)
for j, v in enumerate(syn6): by6[proj_key(v)].append(j)
full_words = []
for i, v in enumerate(syn5):
    key = proj_key((-v) % p)
    for j in by6.get(key, ()):
        # lambda with syn5[i] + lambda * syn6[j] = 0
        nz = np.nonzero(syn6[j])[0][0]; lam = (-int(syn5[i][nz]) * pow(int(syn6[j][nz]), p - 2, p)) % p
        e = (vec([(full_col(5, c), a) for c, a in words[i]]) + lam * vec([(full_col(6, c), a) for c, a in words[j]])) % p
        assert not ((G0b @ e) % p).any() and np.count_nonzero(e) == 8
        full_words.append(("i", e))
paired = sum(1 for i, w in enumerate(words) for j in by6.get(proj_key((-syn5[i]) % p), ()) if words[j] == w)
print(f"(i) cross words: {len(full_words)} (of which (a'|-a') type: {paired}) [{time.time()-T0:.0f}s]", flush=True)

# (ii) structural: two lines with 6 points of Y in total, off-Y part 8 points, full-support dependency of H6
_dirs = sorted({min(tuple(int(v) for v in (c * np.array(d)) % p) for c in range(1, p))
                for d in itertools.product(range(p), repeat=3) if any(d)})
LINES = sorted({tuple(sorted(idx3[tuple(int(v) for v in (P3[a] + c * np.array(d)) % p)] for c in range(p)))
                for a in range(p ** 3) for d in _dirs})
yc = [sum(x in Yset for x in L) for L in LINES]
cand = [li for li, c in enumerate(yc) if c >= 2]
fam = 0; tested = 0; d2words = set()
for a_, b_ in itertools.combinations(cand, 2):
    U = set(LINES[a_]) | set(LINES[b_]); off = sorted(x for x in U if x not in Yset)
    if len(off) != 8: continue
    tested += 1
    ns = np.asarray(_null_space(H6[:, [posY[x] for x in off]], p)) % p
    for v in ns:                                   # with d'' >= 8, any null vector on 8 columns has full support
        if np.count_nonzero(v) == 8:
            nz = int(np.nonzero(v)[0][0]); inv = pow(int(v[nz]), p - 2, p)
            d2words.add(tuple(zip(off, (int(x) * inv % p for x in v))))
print(f"(ii) line-pair 8-sets tested {tested}, d''-words found {len(d2words)} [{time.time()-T0:.0f}s]", flush=True)
exact_ii = None
if len(sys.argv) > 5:                              # exact part (ii) from w8_mitm44v2.py (words as [[point, coeff], ...])
    ex = json.load(open(sys.argv[5])); assert ex["complete"]
    exact_ii = {tuple((int(x), int(a)) for x, a in w) for w in ex["words"]}
    extra = exact_ii - d2words
    print(f"(ii) exact census: {len(exact_ii)} words; line-pair family {len(d2words)} (all contained: {d2words <= exact_ii}); outside it: {len(extra)}", flush=True)
    for w in sorted(extra):
        pts = [tuple(int(v) for v in P3[x]) for x, _ in w]
        print("   extra word points", pts, flush=True)
    d2words = exact_ii
for w in d2words:
    for hyper in (5, 6):
        e = vec([(kpos[idxX[(hyper,) + tuple(int(v) for v in P3[x])]], a) for x, a in w])
        assert not ((G0b @ e) % p).any()
        full_words.append(("ii", e))

# ---- per-output failures (each projective word stands for p-1 vectors)
fails = np.zeros(kk, dtype=np.int64); by_type = defaultdict(int)
for typ, e in full_words:
    hit = (H1 @ e) % p != 0
    fails += hit; by_type[typ] += 1
fails *= (p - 1)
A8_found = (p - 1) * len(full_words)
res = {"t": t, "sample": s, "n": n, "k": kk, "d": 8, "dprime_words": len(words), "cross_words_i": by_type["i"],
       "paired_type": paired, "dpp_words_found_ii_lines": len(d2words), "dpp_words_ii_vectors_in_code": by_type["ii"],
       "A8_found": A8_found, "L8_found": float(fails.mean()), "fails_max": int(fails.max()), "fails_min": int(fails.min()),
       "part_ii_complete": exact_ii is not None, "seconds": time.time() - T0}
print(res, flush=True); json.dump(res, open(out, "w"), indent=1)
