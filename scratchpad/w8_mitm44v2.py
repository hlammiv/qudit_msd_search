"""Exhaustive weight-8 census of ker(H) over F_7, 4+4 meet-in-the-middle, v2: each 4-sum is generated once as a join of
two pair lists bucketed on the first 3 syndrome coordinates (343 passes; pass b holds the left 4-sums whose first 3
coordinates equal b and streams exactly the right 4-sums whose first 3 coordinates equal -b).  Total work = one left
enumeration (C(N,4)(p-1)^3) + one right enumeration (C(N,4)(p-1)^4), unlike v1 which re-streamed everything per pass.
Same invariants as w8_mitm44.py: exact packed syndromes, canonical split (left = 4 smallest columns) recorded once,
hits = C(8,4) x words.  Left 4-sets i<j<k<l with coefficients (1, b, c, d); right with all coefficients.
usage: PYTHONPATH=.:scratchpad python scratchpad/w8_mitm44v2.py selftest | paired SWEEP.json T SAMPLE OUT.json
(paired mode checkpoints to OUT.json.ckpt every 20 passes; rerun the same command to resume.)"""
import sys, json, time, numpy as np
from numba import njit, prange
from w8_mitm44 import _add, _neg, _hash, pack_columns, canonical, brute_weight8, P

SUB = 24
COORDS = int(__import__('os').environ.get('MITM_COORDS', '3'))
if COORDS >= 4: SUB = 21


@njit(cache=True, inline="always")
def _bucket(a0, C):
    r = np.int64(0); m = np.int64(1)
    for t in range(C):
        r += m * np.int64((a0 >> np.uint64(4 * t)) & np.uint64(15)); m *= 7
    return r


@njit(cache=True)
def _pairs(mult, N, normalized, C):
    """Pair list (i<j) with coefficients (1, b) [normalized] or (a, b); returns arrays sorted by (bucket, i)."""
    nc = 6 if normalized else 36
    M = N * (N - 1) // 2 * nc
    s0 = np.empty(M, np.uint64); s1 = np.empty(M, np.uint64); s2 = np.empty(M, np.uint64)
    I = np.empty(M, np.int64); J = np.empty(M, np.int64); CI = np.empty(M, np.int64); CJ = np.empty(M, np.int64); bk = np.empty(M, np.int64)
    q = 0
    for i in range(N):
        for j in range(i + 1, N):
            for ai in range(1, 7):
                if normalized and ai != 1:
                    continue
                for aj in range(1, 7):
                    a0, a1, a2 = _add(mult[ai, i, 0], mult[ai, i, 1], mult[ai, i, 2], mult[aj, j, 0], mult[aj, j, 1], mult[aj, j, 2])
                    s0[q] = a0; s1[q] = a1; s2[q] = a2; I[q] = i; J[q] = j; CI[q] = ai; CJ[q] = aj; bk[q] = _bucket(a0, C)
                    q += 1
    order = np.argsort(bk * N + I, kind="mergesort")
    NB = 7 ** C
    start = np.zeros(NB + 1, np.int64)
    for t in range(M):
        start[bk[t] + 1] += 1
    for t in range(NB):
        start[t + 1] += start[t]
    return s0[order], s1[order], s2[order], I[order], J[order], CI[order], CJ[order], start


@njit(cache=True)
def _diff_table(C):
    NB = 7 ** C
    D = np.zeros((NB, NB), np.int64)
    for x in range(NB):
        for y in range(NB):
            r = 0; m = 1; xx = x; yy = y
            for _ in range(C):
                r += ((xx % 7 - yy % 7) % 7) * m; xx //= 7; yy //= 7; m *= 7
            D[x, y] = r
    return D


@njit(cache=True)
def _neg_bucket(b, C):
    r = 0; m = 1
    for _ in range(C):
        r += ((7 - b % 7) % 7) * m; b //= 7; m *= 7
    return r


@njit(cache=True, parallel=True)
def _left_count(A, G, b, D):
    a0, a1, a2, aI, aJ, aCI, aCJ, aS = A
    g0, g1, g2, gI, gJ, gCI, gCJ, gS = G
    NB = D.shape[0]
    cnt = np.zeros(NB, np.int64)
    for beta in prange(NB):
        gb = D[b, beta]; c = 0
        lo, hi = gS[gb], gS[gb + 1]
        for x in range(aS[beta], aS[beta + 1]):
            jj = aJ[x]
            # G bucket sorted by i: entries with gI > jj form a suffix
            L, R = lo, hi
            while L < R:
                mid = (L + R) // 2
                if gI[mid] > jj: R = mid
                else: L = mid + 1
            c += hi - L
        cnt[beta] = c
    return cnt


@njit(cache=True, parallel=True)
def _left_fill(A, G, b, D, off, keys, pay):
    a0, a1, a2, aI, aJ, aCI, aCJ, aS = A
    g0, g1, g2, gI, gJ, gCI, gCJ, gS = G
    for beta in prange(D.shape[0]):
        gb = D[b, beta]; q = off[beta]
        lo, hi = gS[gb], gS[gb + 1]
        for x in range(aS[beta], aS[beta + 1]):
            jj = aJ[x]
            L, R = lo, hi
            while L < R:
                mid = (L + R) // 2
                if gI[mid] > jj: R = mid
                else: L = mid + 1
            for y in range(L, hi):
                u0, u1, u2 = _add(a0[x], a1[x], a2[x], g0[y], g1[y], g2[y])
                keys[q] = _hash(u0, u1, u2)
                pay[q] = (np.uint64(aI[x]) | (np.uint64(aJ[x]) << np.uint64(9)) | (np.uint64(gI[y]) << np.uint64(18)) |
                          (np.uint64(gJ[y]) << np.uint64(27)) | (np.uint64(aCJ[x]) << np.uint64(36)) |
                          (np.uint64(gCI[y]) << np.uint64(39)) | (np.uint64(gCJ[y]) << np.uint64(42)))
                q += 1


@njit(cache=True)
def _csort(keys, pay):
    """Counting sort by the top SUB bits of the hash; returns (keys, pay, start)."""
    nsub = 1 << SUB; E = keys.shape[0]
    start = np.zeros(nsub + 1, np.int64)
    for t in range(E):
        start[(keys[t] >> np.uint64(64 - SUB)) + np.uint64(1)] += 1
    for t in range(nsub):
        start[t + 1] += start[t]
    pos = start[:-1].copy(); k2 = np.empty_like(keys); p2 = np.empty_like(pay)
    for t in range(E):
        sb = keys[t] >> np.uint64(64 - SUB); q = pos[sb]; pos[sb] += 1
        k2[q] = keys[t]; p2[q] = pay[t]
    return k2, p2, start


@njit(cache=True, parallel=True)
def _right_scan(mult, G, b, D, keys, pay, start, CAP):
    g0, g1, g2, gI, gJ, gCI, gCJ, gS = G
    NB = D.shape[0]; Cc = 0; t_ = NB
    while t_ > 1:
        t_ //= 7; Cc += 1
    nb = _neg_bucket(b, Cc)
    hits = np.zeros(NB, np.int64); nw = np.zeros(NB, np.int64); words = np.zeros((NB, CAP, 16), np.int64)
    for beta in prange(NB):
        qb = D[nb, beta]
        lo, hi = gS[qb], gS[qb + 1]
        for x in range(gS[beta], gS[beta + 1]):
            i = gI[x]; j = gJ[x]
            L, R = lo, hi
            while L < R:
                mid = (L + R) // 2
                if gI[mid] > j: R = mid
                else: L = mid + 1
            for y in range(L, hi):
                k = gI[y]; l = gJ[y]
                u0, u1, u2 = _add(g0[x], g1[x], g2[x], g0[y], g1[y], g2[y])
                v0, v1, v2 = _neg(u0, u1, u2)
                h = _hash(v0, v1, v2); sb = h >> np.uint64(64 - SUB)
                for q in range(start[sb], start[sb + 1]):
                    if keys[q] != h:
                        continue
                    z = pay[q]
                    li = np.int64(z & np.uint64(511)); lj = np.int64((z >> np.uint64(9)) & np.uint64(511))
                    lk = np.int64((z >> np.uint64(18)) & np.uint64(511)); ll = np.int64((z >> np.uint64(27)) & np.uint64(511))
                    if li == i or li == j or li == k or li == l or lj == i or lj == j or lj == k or lj == l:
                        continue
                    if lk == i or lk == j or lk == k or lk == l or ll == i or ll == j or ll == k or ll == l:
                        continue
                    cj = np.int64((z >> np.uint64(36)) & np.uint64(7)); ck = np.int64((z >> np.uint64(39)) & np.uint64(7)); cl = np.int64((z >> np.uint64(42)) & np.uint64(7))
                    w0, w1, w2 = _add(mult[1, li, 0], mult[1, li, 1], mult[1, li, 2], mult[cj, lj, 0], mult[cj, lj, 1], mult[cj, lj, 2])
                    w0, w1, w2 = _add(w0, w1, w2, mult[ck, lk, 0], mult[ck, lk, 1], mult[ck, lk, 2])
                    w0, w1, w2 = _add(w0, w1, w2, mult[cl, ll, 0], mult[cl, ll, 1], mult[cl, ll, 2])
                    if w0 != v0 or w1 != v1 or w2 != v2:
                        continue
                    hits[beta] += 1
                    if ll < i:
                        r_ = nw[beta]
                        if r_ < CAP:
                            words[beta, r_, 0] = li; words[beta, r_, 1] = lj; words[beta, r_, 2] = lk; words[beta, r_, 3] = ll
                            words[beta, r_, 4] = i; words[beta, r_, 5] = j; words[beta, r_, 6] = k; words[beta, r_, 7] = l
                            words[beta, r_, 8] = 1; words[beta, r_, 9] = cj; words[beta, r_, 10] = ck; words[beta, r_, 11] = cl
                            words[beta, r_, 12] = gCI[x]; words[beta, r_, 13] = gCJ[x]; words[beta, r_, 14] = gCI[y]; words[beta, r_, 15] = gCJ[y]
                        nw[beta] += 1
    return hits, nw, words


def census(H, CAP=64, log=print, passes=None, C=None, ckpt=None, ckpt_every=20):
    """``ckpt``: path of a checkpoint (json: fingerprint of (H, C), next pass, words, hits, left); written every
    ``ckpt_every`` passes and resumed from when it matches.  Only for the full pass range (passes=None)."""
    import os, hashlib
    C = COORDS if C is None else C; NB = 7 ** C
    H = np.asarray(H, dtype=np.int64) % P; r, N = H.shape
    assert r >= 3
    mult = pack_columns(H); T0 = time.time()
    A = _pairs(mult, N, True, C); G = _pairs(mult, N, False, C); D = _diff_table(C)
    log(f"pairs: A {A[0].shape[0]}, G {G[0].shape[0]} [{time.time() - T0:.0f}s]")
    allw = []; tot_hits = 0; tot_left = 0; b0 = 0
    fp = hashlib.sha256(np.ascontiguousarray(H).tobytes() + bytes([C])).hexdigest()
    if ckpt and passes is None and os.path.exists(ckpt):
        ck = json.load(open(ckpt))
        if ck.get("fingerprint") == fp:
            allw, tot_hits, tot_left, b0 = ck["words"], ck["hits"], ck["left"], ck["next_pass"]
            log(f"resuming from {ckpt} at pass {b0}/{NB} ({len(allw)} words, {tot_hits} hits)")
        else:
            log(f"ignoring {ckpt}: different matrix or bucketing")

    def save(nxt):
        tmp = ckpt + ".tmp"
        json.dump({"fingerprint": fp, "next_pass": nxt, "words": allw, "hits": tot_hits, "left": tot_left}, open(tmp, "w"))
        os.replace(tmp, ckpt)

    for b in (range(b0, NB) if passes is None else passes):
        t0 = time.time()
        cnt = _left_count(A, G, b, D); off = np.zeros(NB + 1, np.int64); off[1:] = np.cumsum(cnt); E = int(off[-1])
        keys = np.empty(E, np.uint64); pay = np.empty(E, np.uint64)
        _left_fill(A, G, b, D, off, keys, pay)
        keys, pay, start = _csort(keys, pay)
        t1 = time.time()
        hits, nw, words = _right_scan(mult, G, b, D, keys, pay, start, CAP)
        assert (nw <= CAP).all(), "word buffer overflow: raise CAP"
        for beta in np.nonzero(nw)[0]:
            for q in range(nw[beta]):
                allw.append(words[beta, q].tolist())
        tot_hits += int(hits.sum()); tot_left += E
        if b % (20 * 7 ** (C - 3)) == 0 or b == NB - 1:
            log(f"pass {b + 1}/{NB}: left {E} [{t1 - t0:.1f}s] right [{time.time() - t1:.1f}s]; total left {tot_left}, hits {tot_hits}, words {len(allw)} [{time.time() - T0:.0f}s]")
        del keys, pay
        if ckpt and passes is None and ((b + 1) % ckpt_every == 0 or b == NB - 1):
            save(b + 1)
    return allw, tot_hits, tot_left


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "selftest":
        from math import comb
        rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
        cases = [(9, 16, None), (10, 18, None), (11, 20, None), (40, 20, 11)]
        for r, N, low in cases:
            if low:
                Hs = rng.integers(0, P, (low, N)); H = (rng.integers(0, P, (r, low)) @ Hs) % P
            else:
                Hs = H = rng.integers(0, P, (r, N))
            ws, hits, left = census(H, log=lambda m: None)
            got = {canonical(w) for w in ws}; ref = brute_weight8(Hs)
            ok = got == ref and len(got) == len(ws) and hits == 70 * len(got) and left == comb(N, 4) * 216
            print(f"selftest r={r} N={N}: v2 {len(got)} words, brute {len(ref)}, hits {hits} = 70 x {hits // 70}, left {left} = C(N,4)*216: {ok}", flush=True)
            assert ok
    elif mode == "paired":
        from qmsd.reedmuller import rm_generator
        from qmsd.triorthogonal import shorten_matrix
        from qmsd.structured_ad import _rref
        from math import comb
        src, t, s, out = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
        Y = [d for d in json.load(open(src))["draws"] if d["t"] == t and d["sample"] == s][0]["Y"]; Yset = set(Y)
        offY = [x for x in range(P ** 3) if x not in Yset]
        R6 = np.asarray(rm_generator(6, 3, P), dtype=np.int64) % P
        R, pv = _rref(np.asarray(shorten_matrix(R6, Y, P), dtype=np.int64) % P, P); H = R[:len(pv)]
        print(f"t={t}: H {H.shape}", flush=True)
        ws, hits, left = census(H, log=lambda m: print(m, flush=True), ckpt=out + ".ckpt")
        got = sorted({canonical(w) for w in ws})
        N = H.shape[1]
        assert len(got) == len(ws) and hits == 70 * len(got) and left == comb(N, 4) * 216, (len(got), len(ws), hits, left)
        words = [[(offY[c], a) for c, a in w] for w in got]
        json.dump({"t": t, "sample": s, "complete": True, "method": "4+4 MITM v2 (bucketed pair joins), exact packed syndromes",
                   "n_words": len(words), "hits": hits, "left_entries": left, "words": words}, open(out, "w"))
        print(f"DONE t={t}: {len(words)} projective weight-8 d''-words (complete; left entries {left} = C(N,4)*216)", flush=True)
