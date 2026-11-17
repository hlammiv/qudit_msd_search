"""Exhaustive weight-8 census of ker(H) over F_7 by a 4+4 meet-in-the-middle with exact packed syndromes.

H: r x N (r <= 48, N < 512).  Syndromes are packed 4 bits per coordinate in three uint64 (mod-7 SWAR add), so
equality of packed syndromes is exact (no false positives survive the final compare).  Left entries are 4-subsets
i<j<k<l with coefficients (1, b, c, d) (projective normalisation); right entries are 4-subsets with all (p-1)^4
coefficients; a hit is left + right = 0 with disjoint supports.  A weight-8 vector is recorded only from its
canonical split (left = its 4 smallest columns), so each projective weight-8 vector is recorded exactly once; all
hits are also counted (each projective vector appears C(8,4) = 70 times; a word with fewer than 8 nonzeros never
appears since every coefficient is nonzero).  Memory: the left table is split into B hash passes.
usage: PYTHONPATH=. python scratchpad/w8_mitm44.py MODE ... (see main)"""
import sys, json, time, math, numpy as np, numba as nb
from numba import njit, prange

P = 7
ONES = np.uint64(0x1111111111111111)
SEVENS = np.uint64(0x7777777777777777)
M64 = np.uint64(0xFFFFFFFFFFFFFFFF)


@njit(cache=True, inline="always")
def _red(s):
    m = ((s + np.uint64(0x1111111111111111)) >> np.uint64(3)) & np.uint64(0x1111111111111111)
    return s - np.uint64(7) * m


@njit(cache=True, inline="always")
def _add(a0, a1, a2, b0, b1, b2):
    return _red(a0 + b0), _red(a1 + b1), _red(a2 + b2)


@njit(cache=True, inline="always")
def _neg(a0, a1, a2):
    return _red(np.uint64(0x7777777777777777) - a0), _red(np.uint64(0x7777777777777777) - a1), _red(np.uint64(0x7777777777777777) - a2)


@njit(cache=True, inline="always")
def _hash(a0, a1, a2):
    h = a0 * np.uint64(0x9E3779B97F4A7C15)
    h ^= (h >> np.uint64(29)); h += a1 * np.uint64(0xBF58476D1CE4E5B9)
    h ^= (h >> np.uint64(31)); h += a2 * np.uint64(0x94D049BB133111EB)
    h ^= (h >> np.uint64(32)); h *= np.uint64(0xD6E8FEB86659FD93); h ^= (h >> np.uint64(32))
    return h


def pack_columns(H):
    """mult[a, c] = packed (a * H[:, c]) mod 7, a = 0..6 -> uint64 (7, N, 3)."""
    r, N = H.shape
    assert r <= 48 and N < 512
    out = np.zeros((P, N, 3), dtype=np.uint64)
    for a in range(P):
        V = (a * H) % P
        for c in range(N):
            for i in range(r):
                out[a, c, i // 16] |= np.uint64(int(V[i, c])) << np.uint64(4 * (i % 16))
    return out


@njit(cache=True)
def _left_count(mult, N, B, b, SUBBITS):
    """Per-sub-bucket counts of left entries whose hash falls in pass bucket b (per-i rows, summed later)."""
    nsub = 1 << SUBBITS
    cnt = np.zeros((N, nsub), dtype=np.int32) if False else np.zeros(nsub, dtype=np.int64)
    for i in range(N):
        for j in range(i + 1, N):
            for bj in range(1, 7):
                s0, s1, s2 = _add(mult[1, i, 0], mult[1, i, 1], mult[1, i, 2], mult[bj, j, 0], mult[bj, j, 1], mult[bj, j, 2])
                for k in range(j + 1, N):
                    for bk in range(1, 7):
                        t0, t1, t2 = _add(s0, s1, s2, mult[bk, k, 0], mult[bk, k, 1], mult[bk, k, 2])
                        for l in range(k + 1, N):
                            for bl in range(1, 7):
                                u0, u1, u2 = _add(t0, t1, t2, mult[bl, l, 0], mult[bl, l, 1], mult[bl, l, 2])
                                h = _hash(u0, u1, u2)
                                if (h % np.uint64(B)) == np.uint64(b):
                                    cnt[(h >> np.uint64(64 - SUBBITS))] += 1
    return cnt


@njit(cache=True)
def _left_fill(mult, N, B, b, SUBBITS, start, keys, pay):
    pos = start.copy()
    for i in range(N):
        for j in range(i + 1, N):
            for bj in range(1, 7):
                s0, s1, s2 = _add(mult[1, i, 0], mult[1, i, 1], mult[1, i, 2], mult[bj, j, 0], mult[bj, j, 1], mult[bj, j, 2])
                for k in range(j + 1, N):
                    for bk in range(1, 7):
                        t0, t1, t2 = _add(s0, s1, s2, mult[bk, k, 0], mult[bk, k, 1], mult[bk, k, 2])
                        for l in range(k + 1, N):
                            for bl in range(1, 7):
                                u0, u1, u2 = _add(t0, t1, t2, mult[bl, l, 0], mult[bl, l, 1], mult[bl, l, 2])
                                h = _hash(u0, u1, u2)
                                if (h % np.uint64(B)) == np.uint64(b):
                                    sb = h >> np.uint64(64 - SUBBITS)
                                    q = pos[sb]; pos[sb] += 1
                                    keys[q] = h
                                    pay[q] = (np.uint64(i) | (np.uint64(j) << np.uint64(9)) | (np.uint64(k) << np.uint64(18)) |
                                              (np.uint64(l) << np.uint64(27)) | (np.uint64(bj) << np.uint64(36)) |
                                              (np.uint64(bk) << np.uint64(39)) | (np.uint64(bl) << np.uint64(42)))


@njit(cache=True, parallel=True)
def _right_scan(mult, N, B, b, SUBBITS, start, keys, pay, CAP):
    """Thread-parallel right stream over leading column i.  Returns (hits per i, canonical words per i)."""
    hits = np.zeros(N, dtype=np.int64)
    nw = np.zeros(N, dtype=np.int64)
    words = np.zeros((N, CAP, 16), dtype=np.int64)          # 8 cols + 8 coeffs (left first coeff = 1)
    for i in prange(N):
        for ai in range(1, 7):
            for j in range(i + 1, N):
                for aj in range(1, 7):
                    s0, s1, s2 = _add(mult[ai, i, 0], mult[ai, i, 1], mult[ai, i, 2], mult[aj, j, 0], mult[aj, j, 1], mult[aj, j, 2])
                    for k in range(j + 1, N):
                        for ak in range(1, 7):
                            t0, t1, t2 = _add(s0, s1, s2, mult[ak, k, 0], mult[ak, k, 1], mult[ak, k, 2])
                            for l in range(k + 1, N):
                                for al in range(1, 7):
                                    u0, u1, u2 = _add(t0, t1, t2, mult[al, l, 0], mult[al, l, 1], mult[al, l, 2])
                                    v0, v1, v2 = _neg(u0, u1, u2)
                                    h = _hash(v0, v1, v2)
                                    if (h % np.uint64(B)) != np.uint64(b):
                                        continue
                                    sb = h >> np.uint64(64 - SUBBITS)
                                    for q in range(start[sb], start[sb + 1]):
                                        if keys[q] != h:
                                            continue
                                        x = pay[q]
                                        li = np.int64(x & np.uint64(511)); lj = np.int64((x >> np.uint64(9)) & np.uint64(511))
                                        lk = np.int64((x >> np.uint64(18)) & np.uint64(511)); ll = np.int64((x >> np.uint64(27)) & np.uint64(511))
                                        if li == i or li == j or li == k or li == l or lj == i or lj == j or lj == k or lj == l:
                                            continue
                                        if lk == i or lk == j or lk == k or lk == l or ll == i or ll == j or ll == k or ll == l:
                                            continue
                                        cj = np.int64((x >> np.uint64(36)) & np.uint64(7)); ck = np.int64((x >> np.uint64(39)) & np.uint64(7))
                                        cl = np.int64((x >> np.uint64(42)) & np.uint64(7))
                                        w0, w1, w2 = _add(mult[1, li, 0], mult[1, li, 1], mult[1, li, 2], mult[cj, lj, 0], mult[cj, lj, 1], mult[cj, lj, 2])
                                        w0, w1, w2 = _add(w0, w1, w2, mult[ck, lk, 0], mult[ck, lk, 1], mult[ck, lk, 2])
                                        w0, w1, w2 = _add(w0, w1, w2, mult[cl, ll, 0], mult[cl, ll, 1], mult[cl, ll, 2])
                                        if w0 != v0 or w1 != v1 or w2 != v2:
                                            continue
                                        hits[i] += 1
                                        if ll < i and nw[i] < CAP:          # canonical split: left = 4 smallest columns
                                            r_ = nw[i]
                                            words[i, r_, 0] = li; words[i, r_, 1] = lj; words[i, r_, 2] = lk; words[i, r_, 3] = ll
                                            words[i, r_, 4] = i; words[i, r_, 5] = j; words[i, r_, 6] = k; words[i, r_, 7] = l
                                            words[i, r_, 8] = 1; words[i, r_, 9] = cj; words[i, r_, 10] = ck; words[i, r_, 11] = cl
                                            words[i, r_, 12] = ai; words[i, r_, 13] = aj; words[i, r_, 14] = ak; words[i, r_, 15] = al
                                            nw[i] += 1
                                        elif ll < i:
                                            nw[i] += 1                     # overflow: counted, not stored
    return hits, nw, words


def census(H, B, SUBBITS=24, CAP=4096, log=print):
    H = np.asarray(H, dtype=np.int64) % P; r, N = H.shape
    mult = pack_columns(H); T0 = time.time(); allw = []; tot_hits = 0; tot_words = 0
    for b in range(B):
        t0 = time.time()
        cnt = _left_count(mult, N, B, b, SUBBITS)
        start = np.zeros(len(cnt) + 1, dtype=np.int64); start[1:] = np.cumsum(cnt)
        E = int(start[-1]); keys = np.empty(E, dtype=np.uint64); pay = np.empty(E, dtype=np.uint64)
        _left_fill(mult, N, B, b, SUBBITS, start, keys, pay)
        t1 = time.time()
        hits, nw, words = _right_scan(mult, N, B, b, SUBBITS, start, keys, pay, CAP)
        assert (nw <= CAP).all(), "word buffer overflow: raise CAP"
        for i in range(N):
            for q in range(nw[i]):
                allw.append(words[i, q].tolist())
        tot_hits += int(hits.sum()); tot_words += int(nw.sum())
        log(f"pass {b + 1}/{B}: left {E} entries [{t1 - t0:.0f}s], right [{time.time() - t1:.0f}s], hits {int(hits.sum())}, words {int(nw.sum())}; total words {tot_words} [{time.time() - T0:.0f}s]")
        del keys, pay
    return allw, tot_hits


def canonical(w):
    """(cols sorted, coeffs normalised so the first is 1)."""
    pairs = sorted(zip(w[:8], w[8:]))
    inv = pow(int(pairs[0][1]), P - 2, P)
    return tuple((int(c), int(a) * inv % P) for c, a in pairs)


def brute_weight8(H):
    """Projective weight-8 vectors of ker H by brute force over 8-subsets (small N only)."""
    import itertools
    from qmsd.structured_ad import _null_space
    H = np.asarray(H, dtype=np.int64) % P; N = H.shape[1]; out = set()
    for S in itertools.combinations(range(N), 8):
        ns = np.asarray(_null_space(H[:, list(S)], P)) % P
        if ns.shape[0] == 0 or not (ns != 0).any(axis=0).all(): continue
        for coef in itertools.product(range(P), repeat=ns.shape[0]):
            v = (np.array(coef) @ ns) % P
            if (v != 0).all():
                out.add(canonical(list(S) + v.tolist()))
    return out


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "selftest":                          # random small H vs brute force
        rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
        for (r, N) in ((9, 16), (10, 18), (11, 20)):
            H = rng.integers(0, P, (r, N))
            ws, hits = census(H, 3, SUBBITS=10)
            got = {canonical(w) for w in ws}
            assert len(got) == len(ws), "a word was recorded twice"
            ref = brute_weight8(H)
            for w in got:                           # exact kernel check
                v = np.zeros(N, dtype=np.int64)
                for c, a in w: v[c] = a
                assert not ((H @ v) % P).any()
            print(f"selftest r={r} N={N}: mitm {len(got)} words ({hits} hits, 70 x words = {70 * len(got)}), brute {len(ref)}, equal {got == ref}", flush=True)
            assert got == ref and hits == 70 * len(got)
    elif mode == "paired":                          # paired SWEEP.json T SAMPLE B OUT.json
        from qmsd.reedmuller import points, rm_generator
        from qmsd.triorthogonal import shorten_matrix
        from qmsd.structured_ad import _rref
        src, t, s, B, out = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
        Y = [d for d in json.load(open(src))["draws"] if d["t"] == t and d["sample"] == s][0]["Y"]; Yset = set(Y)
        offY = [x for x in range(P ** 3) if x not in Yset]
        R6 = np.asarray(rm_generator(6, 3, P), dtype=np.int64) % P
        R, pv = _rref(np.asarray(shorten_matrix(R6, Y, P), dtype=np.int64) % P, P); H = R[:len(pv)]
        print(f"t={t}: H {H.shape}, B={B}", flush=True)
        ws, hits = census(H, B, log=lambda m: print(m, flush=True))
        got = sorted({canonical(w) for w in ws})
        assert len(got) == len(ws) and hits == 70 * len(got), (len(got), len(ws), hits)
        words = [[(offY[c], a) for c, a in w] for w in got]
        json.dump({"t": t, "sample": s, "complete": True, "method": "4+4 MITM, exact packed syndromes", "n_words": len(words),
                   "hits": hits, "words": words}, open(out, "w"))
        print(f"DONE t={t}: {len(words)} projective weight-8 d''-words (complete)", flush=True)
