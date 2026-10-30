"""Full distance certificate for a sweep survivor of the p = 7, m = 4, r = 7 sphere-charge code.

A survivor (runs/weighted/plane_sweep/plane_sweep_k*.jsonl, d_plane >= 7, no class-block word of
weight <= 5) is certified in five resumable stages, each appended to OUT.json:

  1. setup        G0 (shortened generator), class block Gc (class columns, zero rows dropped), ranks
  2. block6       no class-block word of weight <= 5 (recheck) and no weight-6 block word (3+3 MITM)
  3. full5        no word of weight <= 5 on the full code (serial MITM d <= 4, parallel 2+3 weight-5)
  4. mixed6       no weight-6 word touching a non-class column: the kernel rows force the non-class
                  part into RM(20,4) (d = 5, full-support weight-5/6 words are line-supported --
                  lemma verified once in scratchpad/mixed_w6_check.py), so enumerate the line words
                  of weight 5/6 inside the non-class points and test extension by <= 1 class column
  5. witness7     an explicit weight-7 word: taken from flat_witness_k*.jsonl (planar) and verified,
                  or found by the 3+4 scan on the class block (early exit; a full scan with no hit
                  leaves d >= 8 open, pending the weight-7 mixed-word lemma)

Stages 2-4 together prove d >= 7; stage 5 makes d = 7 exact.
usage: certify_candidate.py K SEED OUT.json [--skip-full5]

--skip-full5: omit stage 3.  It is implied by stages 2 and 4: a word of weight <= 6 either lies in the
class block (stage 2 excludes weight <= 6 there) or has a non-class part, which the 35 kernel rows
force into RM(20,4) (d = 5), whose weight-5/6 words are line-supported (lemma), and stage 4
enumerates exactly those line words with weight 5 (no class part) or 6 / 5+1 and tests extension.
Stage 3 is therefore an independent cross-check (10 h at 8 throttled threads), not a logical need.
"""
import sys, os, json, time, itertools, glob
import numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.weighted import weighted_support, quadric, class_points
from qmsd.structured_ad import _rank, _null_space
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import (weight_d_exists_fast, weight5_exists_parallel, weight6_witness_3x3,
                             weight7_witness_3x4)

k, seed, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]; SKIP_FULL5 = "--skip-full5" in sys.argv
p, m, r = 7, 4, 7
T0 = time.time()
def log(msg):
    print(f"[{time.time()-T0:7.0f}s] {msg}", flush=True)
state = json.load(open(out)) if os.path.exists(out) else {"k": k, "seed": seed, "stages": {}}
def save():
    json.dump(state, open(out, "w"), indent=1)
def done(stage):
    return stage in state["stages"]

# ---- the draw
rec = None
for f in glob.glob("runs/weighted/plane_sweep/plane_sweep_k*.jsonl"):
    for l in open(f):
        if not l.strip(): continue
        rr = json.loads(l)
        if rr["k"] == k and rr["seed"] == seed and "cols_in_X" in rr: rec = rr
assert rec is not None, "survivor not found in the sweep files"
cols = rec["cols_in_X"]
P = np.asarray(points(m, p), dtype=np.int64)
index = {tuple(int(v) for v in x): i for i, x in enumerate(P)}
g = quadric(p, m, (1, 1, 1, 1)); GX, X, mu = weighted_support(p, m, g, r)
cand = class_points(mu, p, 1); cset = set(cand)
G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p
keep = [l for l in range(GX.shape[1]) if l not in set(cols)]
n = G0.shape[1]
is_class = np.array([keep[j] in cset for j in range(n)])
C_cols = np.nonzero(is_class)[0]; N_cols = np.nonzero(~is_class)[0]
Gc = np.ascontiguousarray(G0[:, C_cols]); Gc = Gc[[t for t in range(Gc.shape[0]) if Gc[t].any()]]
if not done("setup"):
    state.update({"n": n, "cols_in_X": cols, "d_plane": rec["d_plane"], "G0_rows": int(G0.shape[0]), "class_cols": int(len(C_cols)),
                  "nonclass_cols": int(len(N_cols)), "block_rows": int(Gc.shape[0]), "block_rank": int(_rank(Gc, p))})
    state["stages"]["setup"] = {"seconds": round(time.time() - T0)}; save()
log(f"k={k} seed={seed}: n={n}, class block {Gc.shape} rank {state['block_rank']}, d_plane {rec['d_plane']}")

# ---- stage 2: class block, weight <= 6
if not done("block6"):
    t = time.time(); low = None
    for d in (2, 3, 4, 5):
        if weight_d_exists_fast(Gc, p, d): low = d; break
    log(f"block weight<=5: {low}")
    w6 = None
    if low is None:
        w6 = weight6_witness_3x3(Gc, p)
        log(f"block weight-6 witness (3+3): {w6}")
    state["stages"]["block6"] = {"block_w_le5": low, "block_w6_witness_class_idx": w6, "seconds": round(time.time() - t)}; save()
    if low is not None or w6 is not None:
        state["verdict"] = f"d = {low or 6} (class block)"; save(); log(state["verdict"]); sys.exit(0)

# ---- stage 3: full code, weight <= 5 (cross-check; implied by stages 2 + 4, see docstring)
if SKIP_FULL5 and not done("full5"):
    state["stages"]["full5"] = {"skipped": True, "reason": "implied by block6 + mixed6 (see docstring)"}; save(); log("full5 skipped (implied)")
if not done("full5"):
    t = time.time(); low = None
    for d in (2, 3, 4):
        if weight_d_exists_fast(G0, p, d): low = d; break
    log(f"full weight<=4: {low}")
    w5 = None
    if low is None:
        w5 = bool(weight5_exists_parallel(G0, p)); log(f"full weight-5 exists: {w5}")
    state["stages"]["full5"] = {"full_w_le4": low, "full_w5": w5, "seconds": round(time.time() - t)}; save()
    if low is not None or w5:
        state["verdict"] = f"d = {low or 5} (full code)"; save(); log(state["verdict"]); sys.exit(0)

# ---- stage 4: mixed weight-6 words (non-class part a line word of RM(20,4) of weight 5 or 6)
if not done("mixed6"):
    t = time.time()
    R3 = np.asarray(rm_generator(3, m, p), dtype=np.int64) % p
    gpt = [int(X[l]) for l in keep]
    K_dim = G0.shape[0] - _rank(G0[:, C_cols], p)
    w = ((g - 1) * (g - 6)) % p
    Nset = set(gpt[j] for j in N_cols); col_of = {gpt[j]: j for j in range(n)}
    winv = {x: pow(int(w[x]), -1, p) for x in Nset}
    def norm(v):
        v = np.asarray(v) % p; nz = np.nonzero(v)[0]
        return None if len(nz) == 0 else tuple(((v * pow(int(v[nz[0]]), -1, p)) % p).tolist())
    cdir = {}
    for j in C_cols:
        d_ = norm(G0[:, j]); assert d_ is not None and d_ not in cdir; cdir[d_] = j
    dirs = [v for v in itertools.product(range(p), repeat=m) if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1]
    seen = set(); n_lines = 0; cands = 0; hits = []
    for a in range(p ** m):
        for v in dirs:
            pts = tuple(sorted(index[tuple((int(P[a][i]) + s * v[i]) % p for i in range(m))] for s in range(p)))
            if pts in seen: continue
            seen.add(pts); n_lines += 1
            if sum(1 for x in pts if x in Nset) < 5: continue
            ns = np.asarray(_null_space(R3[:, list(pts)], p)) % p
            for coef in itertools.product(range(p), repeat=ns.shape[0]):
                if not any(coef): continue
                u = (np.array(coef) @ ns) % p
                supp = [pts[i] for i in np.nonzero(u)[0]]
                if len(supp) not in (5, 6) or any(x not in Nset for x in supp): continue
                cands += 1
                eN = [(int(u[pts.index(x)]) * winv[x]) % p for x in supp]
                scols = [col_of[x] for x in supp]
                sv = (G0[:, scols] @ np.array(eN)) % p
                if len(supp) == 6:
                    if not sv.any(): hits.append(["weight6_allN", scols, eN])
                else:
                    d_ = norm((-sv) % p)
                    if d_ is None: hits.append(["weight5_N_only", scols, eN])
                    elif d_ in cdir: hits.append(["weight5_N+1class", scols, eN, int(cdir[d_])])
    log(f"mixed weight-6: kernel dim {K_dim} (expect 35), lines {n_lines}, candidates {cands}, hits {len(hits)}")
    state["stages"]["mixed6"] = {"kernel_dim": int(K_dim), "lines": n_lines, "line_word_candidates": cands, "hits": hits, "seconds": round(time.time() - t)}; save()
    if hits:
        state["verdict"] = "d = 6 (mixed word)"; save(); log(state["verdict"]); sys.exit(0)
state["verdict"] = "d >= 7 (block6 + mixed6" + ("" if SKIP_FULL5 else " + full5") + ")"; save(); log(state["verdict"])

# ---- stage 5: weight-7 witness
if not done("witness7"):
    t = time.time(); wit = None; src = None
    for f in glob.glob("runs/weighted/plane_sweep/flat_witness_k*.jsonl"):
        for l in open(f):
            if not l.strip(): continue
            rr = json.loads(l)
            if rr["k"] == k and rr["seed"] == seed and rr.get("witness_weight") == 7:
                wit = rr["witness_G0_cols"]; coef = rr["witness_coeffs"]; src = rr.get("witness_from", "flat")
    if wit is not None:
        vec = np.zeros(n, dtype=np.int64); vec[wit] = coef
        assert ((G0 @ vec) % p == 0).all() and len(wit) == 7
        log(f"weight-7 witness from flat_witness ({src}) verified")
    else:
        log("no stored witness: running the 3+4 scan on the class block (early exit)")
        wc = weight7_witness_3x4(Gc, p, block=8, log=lambda lo, hi, s: log(f"  block cols [{lo},{hi}) {s:.0f}s"))
        if wc is not None:
            wit = [int(C_cols[j]) for j in wc]
            ns = np.asarray(_null_space(G0[:, wit], p)) % p; coef = [int(x) for x in ns[0]]
            vec = np.zeros(n, dtype=np.int64); vec[wit] = coef; assert ((G0 @ vec) % p == 0).all()
            src = "block 3+4 scan"
    state["stages"]["witness7"] = {"witness_G0_cols": wit, "witness_coeffs": coef if wit is not None else None, "source": src, "seconds": round(time.time() - t)}
    if wit is not None:
        state["verdict"] = "d = 7 exact"; state["gamma"] = float(np.log(n / k) / np.log(7))
    else:
        state["verdict"] = "d >= 7; no weight-7 block word (d = 8 needs the weight-7 mixed-word lemma)"
    save(); log(state["verdict"])
