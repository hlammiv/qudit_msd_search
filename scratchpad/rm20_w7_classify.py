"""Weight-7 words of RM_7(20,4) = RM(3,4)^perp: which 7-point sets of F_7^4 carry a FULL-support
annihilator of the cubics?  (Extends step (2) of mixed_w6_check.py from weight 5/6 to weight 7.)

Up to affine equivalence every 7-set of affine span j contains a frame 0, e_1..e_j and 7-(j+1)
further points of the coordinate j-flat, so the sets with span 2, 3, 4 are covered by
C(46,4) + C(339,3) + C(2396,2) = 163 185 + 6 435 689 + 2 869 210 frame sets (span 1 = a full
line).  Two independent methods:

  A. brute force: batched F_7 rank of the 35 x 7 evaluation matrix of RM(3,4) at every set
     (qmsd.pointsets.batched_rank_modp); dependent sets get their null space, a full-support
     test and an affine-invariant classification (collinearity profile);
  B. residual hashing: after row-reducing the frame columns to unit vectors, the extras'
     residuals (rows j+1 .. 34) decide dependence; for span 4 a dependent pair means equal
     residual directions, for span 3 a dependent triple means c' in span(a', b') (6 directions
     per pair), for span 2 a batched rank of the 7 x 4 residual blocks.

Writes runs/weighted/rm20_w7_classification.json (resumable per span).  ~10 min on one core,
< 1 GB.  usage: PYTHONPATH=. python scratchpad/rm20_w7_classify.py [--fast]  (--fast: method B only)
"""
import sys, os, json, time, itertools
from math import comb
import numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.pointsets import (batched_rank_modp, full_support_null_vector, classify_pointset,
                            frame_and_flat, inverse_table)
from qmsd.structured_ad import _null_space

p, m = 7, 4
OUT = "runs/weighted/rm20_w7_classification.json"
FAST = "--fast" in sys.argv
BATCH = 50_000
T0 = time.time()
def log(s): print(f"[{time.time()-T0:6.0f}s] {s}", flush=True)
state = json.load(open(OUT)) if os.path.exists(OUT) else {"p": p, "m": m, "spans": {}}
def save(): json.dump(state, open(OUT, "w"), indent=1)

P = np.asarray(points(m, p), dtype=np.int64)
R3 = np.asarray(rm_generator(3, m, p), dtype=np.int64) % p          # 35 x 2401
assert R3.shape == (35, p ** m)
inv = inverse_table(p)

def combos(n, r):
    it = itertools.combinations(range(n), r)
    while True:
        chunk = np.fromiter(itertools.chain.from_iterable(itertools.islice(it, BATCH)), dtype=np.int64)
        if chunk.size == 0: return
        yield chunk.reshape(-1, r)

def analyse(T):
    """null space + full-support word + label of a dependent 7-set (global point indices)."""
    u = full_support_null_vector(R3[:, T], p)
    ns = np.asarray(_null_space(R3[:, T], p)) % p
    return {"points": [tuple(int(v) for v in P[x]) for x in T], "nullity": int(ns.shape[0]),
            "label": classify_pointset(P[list(T)], p),
            "full_support_word": None if u is None else [int(v) for v in u]}

def frame_reduce(frame):
    """Row-reduce R3 so that the frame columns become e_0..e_j; returns (M, pivot rows)."""
    M = R3.copy(); rows = []
    for i, c in enumerate(frame):
        r = next(t for t in range(i, 35) if M[t, c] != 0)
        M[[i, r]] = M[[r, i]]
        M[i] = (M[i] * inv[M[i, c]]) % p
        for t in range(35):
            if t != i and M[t, c]: M[t] = (M[t] - M[t, c] * M[i]) % p
    return M

def normalise(V):
    """Rows scaled so that the first non-zero entry is 1 (zero rows stay zero)."""
    V = V % p
    nz = V != 0
    first = np.argmax(nz, axis=1)
    lead = V[np.arange(len(V)), first]
    return (V * inv[lead][:, None]) % p

for j in (2, 3, 4):
    key = f"span{j}"
    if key in state["spans"] and state["spans"][key].get("done"): log(f"{key}: done, skipping"); continue
    frame, flat = frame_and_flat(m, j, p)
    others = np.array([x for x in flat if x not in set(frame)], dtype=np.int64)
    extra = 7 - (j + 1)
    rec = {"frame": [tuple(int(v) for v in P[x]) for x in frame], "flat_points": len(flat), "extras": int(extra)}
    # ---- method B: residual hashing
    t = time.time()
    M = frame_reduce(frame)
    resid = M[j + 1:, :][:, others].T.copy()                         # (n_others, 34 - j)
    coefs = M[:j + 1, :][:, others].T.copy()                         # frame coefficients of each extra
    zero_res = int((~resid.any(axis=1)).sum())
    keys = {}
    for i, row in enumerate(normalise(resid)):
        keys.setdefault(row.tobytes(), []).append(i)
    prop_pairs = sum(len(v) * (len(v) - 1) // 2 for v in keys.values())
    rec["residual"] = {"dim": int(34 - j), "zero_residuals": zero_res, "proportional_pairs": prop_pairs}
    B_dep = []                                                       # dependent sets found by method B (sorted tuples of other-indices)
    if j == 4:
        for v in keys.values():
            for a, b in itertools.combinations(v, 2): B_dep.append((a, b))
    elif j == 3:
        assert zero_res == 0 and prop_pairs == 0
        pos = {kk: v[0] for kk, v in keys.items()}
        found = set()
        for a, b in itertools.combinations(range(len(others)), 2):
            cand = normalise(np.array([(resid[a] + beta * resid[b]) % p for beta in range(1, p)]))
            for row in cand:
                c = pos.get(row.tobytes())
                if c is not None: found.add(tuple(sorted((a, b, c))))
        B_dep = sorted(found)
    else:
        n_o = len(others); dep = []
        for chunk in combos(n_o, extra):
            A = resid[chunk]                                         # (B, 4, 7)
            rk = batched_rank_modp(A, p)
            dep.extend(tuple(int(v) for v in row) for row in chunk[rk < extra])
        B_dep = dep
    rec["residual"]["dependent_sets"] = len(B_dep)
    rec["residual"]["seconds"] = round(time.time() - t, 1)
    log(f"{key}: method B (residual): {len(B_dep)} dependent sets, zero residuals {zero_res}, proportional pairs {prop_pairs} [{rec['residual']['seconds']}s]")
    # ---- method A: brute force rank of the 35 x 7 matrices
    if not FAST:
        t = time.time(); n_sets = 0; A_dep = []
        total = comb(len(others), extra)                                  # 163 185 / 6 435 689 / 2 869 210
        for chunk in combos(len(others), extra):
            T = np.concatenate([np.broadcast_to(np.array(frame), (len(chunk), j + 1)), others[chunk]], axis=1)   # (B, 7)
            A = R3.T[T]                                              # (B, 7, 35)
            rk = batched_rank_modp(A, p)
            n_sets += len(chunk)
            for row in chunk[rk < 7]: A_dep.append(tuple(int(v) for v in row))
            if n_sets % (20 * BATCH) == 0 or n_sets == total:
                log(f"  {key}: {n_sets}/{total} sets, dependent so far {len(A_dep)}")
        rec["brute"] = {"sets": n_sets, "dependent_sets": len(A_dep), "seconds": round(time.time() - t)}
        assert n_sets == total
        assert sorted(A_dep) == sorted(tuple(x) for x in B_dep), "methods A and B disagree"
        log(f"{key}: method A (brute force): {n_sets} sets, {len(A_dep)} dependent; agrees with B [{rec['brute']['seconds']}s]")
    # ---- classify the dependent sets
    t = time.time(); labels = {}; full = []; nullities = {}
    for sset in B_dep:
        T = list(frame) + [int(others[i]) for i in sset]
        info = analyse(T)
        labels[info["label"]] = labels.get(info["label"], 0) + 1
        nullities[info["nullity"]] = nullities.get(info["nullity"], 0) + 1
        if info["full_support_word"] is not None: full.append(info)
    rec["dependent_labels"] = labels; rec["dependent_nullities"] = {str(k): v for k, v in nullities.items()}
    rec["full_support_sets"] = len(full)
    rec["full_support_labels"] = {}
    for info in full: rec["full_support_labels"][info["label"]] = rec["full_support_labels"].get(info["label"], 0) + 1
    rec["full_support_examples"] = full[:20]
    rec["done"] = True
    state["spans"][key] = rec; save()
    log(f"{key}: dependent labels {labels}; full-support sets {len(full)} {rec['full_support_labels']} [{time.time()-t:.0f}s]")

# ---- span 1: the full line (RM(2,1) words), counted exactly
line = [i for i, x in enumerate(P) if not x[1:].any()]
ns = np.asarray(_null_space(R3[:, line], p)) % p
seen = set(); wt = {}
for coef in itertools.product(range(p), repeat=ns.shape[0]):
    if not any(coef): continue
    u = (np.array(coef) @ ns) % p
    f = np.nonzero(u)[0][0]; u = tuple(((u * inv[u[f]]) % p).tolist())
    if u in seen: continue
    seen.add(u); w = sum(1 for v in u if v); wt[w] = wt.get(w, 0) + 1
n_lines = p ** (m - 1) * (p ** m - 1) // (p - 1)
state["span1_line"] = {"nullity": int(ns.shape[0]), "projective_words_by_weight": wt, "lines": n_lines,
                       "codewords_by_weight": {w: c * n_lines * (p - 1) for w, c in wt.items()}}
save()
log(f"line: nullity {ns.shape[0]}, projective words by weight {wt}; {n_lines} lines")
log("done: " + OUT)
