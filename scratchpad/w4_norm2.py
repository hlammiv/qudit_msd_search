"""Exhaustive weight-4 test with a scalar-normalised LEFT table (2+2 split): left pairs carry
coefficient patterns (1, c), right pairs all (c3, c4), so the table is (p-1)x smaller than
mindist_nb's full left table (needed at n ~ 12 000 where the full table does not fit).  Every hash
match is verified exactly.  Also usable as a plain weight-<=4 check: run weight 2, 3 first with
mindist_nb (cheap).
usage: w4_norm2.py CAND.json   (rm_m4_candidate.py output; appends stage "4" to it)"""
import sys, json, time, numpy as np, numba as nb
from qmsd.reedmuller import rm_generator
from qmsd.triorthogonal import shorten_matrix

@nb.njit(cache=True)
def _codes_pairs_norm(H, p, pw):
    r, n = H.shape; N = (n * (n - 1) // 2) * (p - 1)
    codes = np.empty(N, dtype=np.uint64); supp = np.empty((N, 2), dtype=np.int32); coef = np.empty(N, dtype=np.int8)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            for c in range(1, p):
                h = np.uint64(0)
                for t in range(r):
                    v = (H[t, i] + c * H[t, j]) % p
                    h = h + pw[t] * np.uint64(v)
                codes[idx] = h; supp[idx, 0] = i; supp[idx, 1] = j; coef[idx] = c; idx += 1
    return codes, supp, coef

@nb.njit(parallel=True, cache=True)
def _scan_right_pairs(H, p, pw, codes_sorted, order, supp, coef, nthreads):
    r, n = H.shape; N = codes_sorted.shape[0]
    found = np.zeros((nthreads, 8), dtype=np.int64); hits = np.zeros(nthreads, dtype=np.int64)
    for k in nb.prange(n):
        tid = nb.get_thread_id() if hasattr(nb, "get_thread_id") else 0
        for l in range(k + 1, n):
            for c3 in range(1, p):
                for c4 in range(1, p):
                    h = np.uint64(0)
                    for t in range(r):
                        v = (-(c3 * H[t, k] + c4 * H[t, l])) % p
                        h = h + pw[t] * np.uint64(v)
                    # binary search
                    lo = 0; hi = N
                    while lo < hi:
                        mid = (lo + hi) >> 1
                        if codes_sorted[mid] < h: lo = mid + 1
                        else: hi = mid
                    m = lo
                    while m < N and codes_sorted[m] == h:
                        q = order[m]; i = supp[q, 0]; j = supp[q, 1]
                        if i != k and i != l and j != k and j != l:
                            ok = True
                            for t in range(r):
                                if (H[t, i] + coef[q] * H[t, j] + c3 * H[t, k] + c4 * H[t, l]) % p != 0:
                                    ok = False; break
                            if ok and hits[tid] == 0:
                                hits[tid] = 1; found[tid, 0] = i; found[tid, 1] = j; found[tid, 2] = k; found[tid, 3] = l
                                found[tid, 4] = 1; found[tid, 5] = coef[q]; found[tid, 6] = c3; found[tid, 7] = c4
                        m += 1
    return hits, found

def weight4_exists_norm2(H, p, nthreads):
    H = np.ascontiguousarray(np.asarray(H, dtype=np.int64) % p); r = H.shape[0]
    pw = np.array([pow(p, t, 2 ** 64) for t in range(r)], dtype=np.uint64)
    t0 = time.time(); codes, supp, coef = _codes_pairs_norm(H, p, pw)
    order = np.argsort(codes, kind="stable").astype(np.int64); codes_sorted = codes[order]
    print(f"  left table {len(codes):.3e} entries built+sorted in {time.time()-t0:.0f}s", flush=True)
    t1 = time.time(); hits, found = _scan_right_pairs(H, p, pw, codes_sorted, order, supp, coef, nthreads)
    print(f"  right scan {time.time()-t1:.0f}s", flush=True)
    for tid in range(nthreads):
        if hits[tid]:
            return [int(x) for x in found[tid, :4]], [int(x) for x in found[tid, 4:]]
    return None, None

if __name__ == "__main__":
    f = sys.argv[1]; c = json.load(open(f)); p, m, r = c["p"], c["m"], c["r"]
    cols0 = [x - 1 for x in c["puncture_columns_1indexed"]]
    G = np.asarray(rm_generator(r, m, p), dtype=np.int64) % p
    G0 = np.asarray(shorten_matrix(G, cols0, p), dtype=np.int64) % p; G0 = G0[[t for t in range(G0.shape[0]) if G0[t].any()]]
    nthreads = nb.get_num_threads(); print(f"G0 {G0.shape}, {nthreads} threads", flush=True)
    t = time.time(); cols, coefs = weight4_exists_norm2(G0, p, nthreads)
    c.setdefault("stages", {})["4"] = {"exists": cols is not None, "witness_cols": cols, "witness_coeffs": coefs, "seconds": round(time.time() - t), "method": "2+2 normalised left (w4_norm2.py)"}
    if cols is None: c["d_lower"] = 5
    else: c["d"] = 4
    json.dump(c, open(f, "w"), indent=1); print(f"weight 4 exists: {cols is not None} {cols} {coefs} ({time.time()-t:.0f}s)", flush=True)
