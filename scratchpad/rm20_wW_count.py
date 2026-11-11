"""Exact number A_W of weight-W words of RM_7(20,4) = RM(3,4)^perp on F_7^4, W = 8 or 9 (generalises
rm20_w8_count.py; W = 8 must reproduce A_8 = 436 272 744 600).

Frame sets: frame 0, e_1..e_j plus e = W-(j+1) extras of the coordinate j-flat; extras' residuals after the frame
columns are row-reduced decide dependence (rank < e).  A full-support word needs every point in a circuit.

  span 2: all C(46, e) sets, batched rank of the residual blocks, then a batched coloop test.
  span 3, 4: a full-support set S either
     (a) contains a dependent (e-1)-subset Q of extras and its last extra x lies in span(resid[Q]) (x not a coloop),
         found by looking up the directions of span(resid[Q]) in the table of single residual directions; or
     (b) has every (e-1)-subset of extras independent, so its extras form a circuit: a split into a pair and the
         rest gives  a + beta b  ~  (combination of the other e-2 with all coefficients non-zero);
         e = 4: pair-pair collisions; e = 5: pair-triple collisions.
     The dependent (e-1)-subsets come from the same recursion one level down (e-1 = 3 in span 4: pair-single
     collisions; e-1 = 4 in span 3: (a)+(b) at e = 4).
Each candidate is checked exactly: null space of the 35 x W evaluation matrix, number of full-support null vectors.

Orbit count: A_W = sum over frame sets R of (#full-support vectors on R) N_j / F(R) (see rm20_w8_count.py), with the
check sum_R N_2/F(R) = 139 650 C(49, W) over all span-2 frame sets.

usage: PYTHONPATH=. python scratchpad/rm20_wW_count.py W   -> runs/weighted/rm20_w{W}_count_v2.json
"""
import sys, os, json, time, itertools
from math import comb
from fractions import Fraction
import numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.pointsets import batched_rank_modp, frame_and_flat, inverse_table, direction_keys, _lexsort_rows, _find_rows
from qmsd.structured_ad import _null_space

W = int(sys.argv[1]); p, m = 7, 4
OUT = f"runs/weighted/rm20_w{W}_count_v2.json"
T0 = time.time()
def log(s): print(f"[{time.time()-T0:6.0f}s] {s}", flush=True)
state = json.load(open(OUT)) if os.path.exists(OUT) else {"p": p, "m": m, "weight": W, "spans": {}}
def save(): json.dump(state, open(OUT, "w"), indent=1)
P = np.asarray(points(m, p), dtype=np.int64)
R3 = np.asarray(rm_generator(3, m, p), dtype=np.int64) % p
inv = inverse_table(p); NPTS = p ** m

def N_tuples(j):
    v = NPTS
    for i in range(j): v *= NPTS - p ** i
    return v

def frames_in(T, j):
    pts = P[list(T)]
    tup = np.array(list(itertools.permutations(range(len(T)), j + 1)), dtype=np.int64)
    D = (pts[tup[:, 1:]] - pts[tup[:, :1]]) % p
    return int((batched_rank_modp(D, p) == j).sum())

def null_info(T):
    ns = np.asarray(_null_space(R3[:, list(T)], p)) % p
    k = ns.shape[0]
    if k == 0: return 0, 0
    assert k <= 7, k
    coefs = np.array(list(itertools.product(range(p), repeat=k)), dtype=np.int64)
    return k, int((((coefs @ ns) % p) != 0).all(axis=1).sum())

def frame_reduce(frame):
    M = R3.copy()
    for i, c in enumerate(frame):
        r = next(t for t in range(i, 35) if M[t, c] != 0)
        M[[i, r]] = M[[r, i]]; M[i] = (M[i] * inv[M[i, c]]) % p
        for t in range(35):
            if t != i and M[t, c]: M[t] = (M[t] - M[t, c] * M[i]) % p
    return M

def span_dirs(V):
    """all projective directions of the row span of V (k x D), as an array of representatives"""
    V = np.asarray(V) % p; k = V.shape[0]
    coefs = np.array([c for c in itertools.product(range(p), repeat=k) if any(c)], dtype=np.int64)
    S = (coefs @ V) % p
    S = S[S.any(axis=1)]
    return S

class SingleTable:
    def __init__(self, resid):
        self.K = direction_keys(resid, p); self.o = _lexsort_rows(self.K); self.Ks = self.K[self.o]
        assert (self.Ks[1:] != self.Ks[:-1]).any(axis=1).all(), "proportional residual pair"
    def find(self, V):
        pos = _find_rows(self.Ks, direction_keys(V, p))
        return np.where(pos >= 0, self.o[np.maximum(pos, 0)], -1)

def pair_table(resid, chunk_rows=2_000_000):
    """keys of a + beta b (a < b, beta = 1..p-1), with (a, b); built per a so the vectors are never all in memory"""
    n = len(resid); Ks, A, B = [], [], []
    for a in range(n - 1):
        bs = np.arange(a + 1, n)
        for beta in range(1, p):
            V = (resid[a][None, :] + beta * resid[bs]) % p
            Ks.append(direction_keys(V, p)); A.append(np.full(len(bs), a, np.int32)); B.append(bs.astype(np.int32))
    K = np.concatenate(Ks); A = np.concatenate(A); B = np.concatenate(B)
    o = _lexsort_rows(K)
    return K[o], A[o], B[o]

def dependent_subsets(resid, size, table):
    """all `size`-subsets (sorted tuples) of extras with dependent residuals, size = 2, 3, 4 (no zero residuals,
    no proportional pairs assumed)"""
    if size == 2:
        return set()
    lower = dependent_subsets(resid, size - 1, table) if size > 3 else set()
    out = set()
    n = len(resid)
    # (a) contain a dependent (size-1)-subset Q: any extra point completes it to a dependent set
    for Q in lower:
        for x in range(n):
            if x not in Q: out.add(tuple(sorted(Q + (x,))))
    # (b) circuits on all `size` extras
    K, A, B = table["pairs"]
    if size == 3:                                     # a + beta b ~ c
        st = table["single"]
        for a in range(n - 1):
            bs = np.arange(a + 1, n)
            for beta in range(1, p):
                c = st.find((resid[a][None, :] + beta * resid[bs]) % p)
                for i in np.nonzero(c >= 0)[0]:
                    if c[i] not in (a, bs[i]): out.add(tuple(sorted((a, int(bs[i]), int(c[i])))))
    elif size == 4:                                   # a + beta b ~ c + delta d, disjoint pairs: equal keys in K
        same = np.zeros(len(K), bool); same[1:] = (K[1:] == K[:-1]).all(axis=1)
        starts = np.nonzero(~same)[0]; ends = np.append(starts[1:], len(K))
        multi = np.nonzero(ends - starts > 1)[0]
        for g in multi:
            s, e = starts[g], ends[g]
            for i in range(s, e):
                for j2 in range(i + 1, e):
                    q = {int(A[i]), int(B[i]), int(A[j2]), int(B[j2])}
                    if len(q) == 4: out.add(tuple(sorted(q)))
    return out

def circuits5(resid, table):
    """5-subsets of extras forming a circuit: a + beta b ~ c + gamma d + delta e (all non-zero, disjoint)"""
    K, A, B = table["pairs"]; n = len(resid); out = set()
    gd = np.array([(g, d) for g in range(1, p) for d in range(1, p)], dtype=np.int64)
    for c in range(n - 2):
        ii, jj = np.triu_indices(n - c - 1, 1); d_ = ii + c + 1; e_ = jj + c + 1
        for g, dl in gd:
            V = (resid[c][None, :] + g * resid[d_] + dl * resid[e_]) % p
            pos = _find_rows(K, direction_keys(V, p))
            for t in np.nonzero(pos >= 0)[0]:
                # every pair with this key (the table may hold several)
                k0 = pos[t]; key = K[k0]
                lo = k0
                while lo > 0 and (K[lo - 1] == key).all(): lo -= 1
                hi = k0
                while hi + 1 < len(K) and (K[hi + 1] == key).all(): hi += 1
                for r in range(lo, hi + 1):
                    q = {int(A[r]), int(B[r]), c, int(d_[t]), int(e_[t])}
                    if len(q) == 5: out.add(tuple(sorted(q)))
        if c % 50 == 0: log(f"    circuits5: c={c}/{n}, found {len(out)}")
    return out

def completions(resid, base_sets, st):
    """sets Q + {x} with resid[x] in span(resid[Q]) (x not a coloop), for dependent Q"""
    out = set()
    for Q in base_sets:
        S = span_dirs(resid[list(Q)])
        xs = st.find(S)
        for x in set(int(v) for v in xs if v >= 0):
            if x not in Q: out.add(tuple(sorted(Q + (x,))))
    return out

for j in (2, 3, 4):
    key = f"span{j}"
    if key in state["spans"] and state["spans"][key].get("done"):
        log(f"{key}: done"); continue
    t = time.time()
    frame, flat = frame_and_flat(m, j, p)
    others = np.array([x for x in flat if x not in set(frame)], dtype=np.int64)
    e = W - (j + 1)
    M = frame_reduce(frame); resid = M[j + 1:, :][:, others].T.copy() % p
    assert resid.any(axis=1).all()
    rec = {"extras": e, "others": len(others)}
    if j == 2:
        cand = []; it = itertools.combinations(range(len(others)), e); nsets = 0
        while True:
            ch = np.fromiter(itertools.chain.from_iterable(itertools.islice(it, 50_000)), dtype=np.int64)
            if ch.size == 0: break
            ch = ch.reshape(-1, e); nsets += len(ch)
            rk = batched_rank_modp(resid[ch], p)
            dep = ch[rk < e]
            if len(dep) == 0: continue
            # coloop test on the full evaluation matrices: every point removable without rank loss
            T = np.concatenate([np.broadcast_to(np.array(frame), (len(dep), j + 1)), others[dep]], axis=1)
            full_rank = batched_rank_modp(R3.T[T], p)
            ok = np.ones(len(dep), bool)
            for i in range(W):
                sub = np.delete(T, i, axis=1)
                ok &= batched_rank_modp(R3.T[sub], p) == full_rank
            cand.extend(tuple(int(v) for v in row) for row in dep[ok])
        rec["sets"] = nsets
        cands = set(cand)
    else:
        st = SingleTable(resid)
        log(f"{key}: building pair table ({len(others)} extras)")
        table = {"single": st, "pairs": pair_table(resid)}
        log(f"{key}: pair table {len(table['pairs'][0])} rows [{time.time()-t:.0f}s]")
        lower = dependent_subsets(resid, e - 1, table)
        rec["dependent_subsets_e-1"] = len(lower)
        log(f"{key}: {len(lower)} dependent {e-1}-subsets")
        cands = completions(resid, lower, st)
        rec["completions"] = len(cands)
        if e in (3, 4):
            circ = {q for q in dependent_subsets(resid, e, table) if not any(tuple(sorted(set(q) - {x})) in lower for x in q)}
        else:
            circ = circuits5(resid, table)
        rec["circuits"] = len(circ)
        cands |= circ
    log(f"{key}: {len(cands)} candidates [{time.time()-t:.0f}s]")
    Nj = N_tuples(j); total = Fraction(0); nfull = 0; hist = {}; examples = []
    for sset in sorted(cands):
        T = list(frame) + [int(others[i]) for i in sset]
        k, full = null_info(T)
        if full == 0: continue
        F = frames_in(T, j)
        total += Fraction(Nj * full, F); nfull += 1
        lab = f"nullity{k}:full{full}:frames{F}"; hist[lab] = hist.get(lab, 0) + 1
        if len(examples) < 10: examples.append([tuple(int(v) for v in P[x]) for x in T])
    assert total.denominator == 1
    rec.update({"candidates": len(cands), "full_support_sets": nfull, "types": hist, "examples": examples,
                "words": int(total), "N_j": Nj, "seconds": round(time.time() - t), "done": True})
    state["spans"][key] = rec; save()
    log(f"{key}: full-support sets {nfull}, types {hist}, words {int(total)}")

if "orbit_check_span2" not in state:
    frame, flat = frame_and_flat(m, 2, p); others = [x for x in flat if x not in set(frame)]
    col = set()
    for tr in itertools.combinations(sorted(flat), 3):
        D = (P[list(tr[1:])] - P[tr[0]]) % p
        if batched_rank_modp(D[None], p)[0] < 2: col.add(tr)
    acc = Fraction(0); Nj = N_tuples(2); nT = W * (W - 1) * (W - 2)
    for ext in itertools.combinations(others, W - 3):
        T = sorted(list(frame) + list(ext))
        nc = sum(1 for tr in itertools.combinations(T, 3) if tr in col)
        acc += Fraction(Nj, nT - 6 * nc)
    state["orbit_check_span2"] = {"ok": acc == 139650 * comb(49, W), "sum": str(acc), "expected": 139650 * comb(49, W)}
    save(); log(f"orbit check: {state['orbit_check_span2']['ok']}")
state[f"A{W}_RM20"] = sum(v["words"] for v in state["spans"].values()); save()
log(f"A_{W}(RM_7(20,4)) = {state[f'A{W}_RM20']}")
