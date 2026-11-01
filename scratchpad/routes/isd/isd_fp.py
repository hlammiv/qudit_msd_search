"""Information-set decoding (ISD) over F_p as a randomised low-weight-word finder for ker(H).

Setting.  H = G0 is an r x n matrix over F_p with r = rank(H) tiny (15-60) and n large
(400-12 000).  We want e in F_p^n with H e = 0 and wt(e) <= w.  H plays the role of the
parity-check matrix of the code C = ker(H), which has redundancy r and rate 1 - r/n > 0.99.

The regime matters.  Classical ISD (Prange / Lee-Brickell / Stern / Dumer) picks a random
information set of C, i.e. a random set B of r independent columns of H, eliminates
H -> [Q | I_r] (columns A = complement of B first), and uses the fact that for ANY choice of e_A
the redundancy part e_B = -Q e_A is free: a word is found iff its support puts only a small
weight t on A and the rest (w - t) on B.  ISD beats brute force because a random r-set catches a
constant fraction of the weight when r is a constant fraction of n.  Here r/n ~ 0.002-0.05, so a
weight-w word has on average w r/n ~ 0.01-0.3 of its support in B: with probability
(1 - r/n)^w ~ 0.75-0.99 the information set absorbs NOTHING and ISD degenerates to "enumerate
weight-w patterns on A and check Q e_A = 0", i.e. to the exhaustive meet-in-the-middle (MITM)
scan the repo already runs (qmsd/mindist_nb.py).  The per-iteration success probabilities below
quantify this: every ISD variant is either (a) an exhaustive MITM in disguise (then it is a
certificate, but costs C(n, w/2)(p-1)^{w/2-1}) or (b) a sub-exhaustive sampler whose per-word
success probability is a product of (r/n) or (z/n) factors, so small that "not found after N
iterations" says nothing useful about d once the exhaustive scan is out of reach.

What is implemented (all return an explicitly verified word or None):
  * ``mitm_exhaustive``: the projective-key MITM (both halves with (p-1)^{t-1} normalised
    patterns, keyed by the syndrome up to scaling) on a column set.  Exhaustive = a proof of
    absence on that column set.  Memory/time ~ C(s, ceil(w/2)) (p-1)^{ceil(w/2)-1}.
  * ``stern_iteration``: p-ary Stern/Dumer: random information set B, random disjoint windows
    Z1, Z2 in A of sizes z1, z2, random set of ell collision rows; lists of weight-t1 / weight-t2
    patterns on Z1 / Z2 keyed by the projective ell-part of Q c; matches are completed on the
    remaining r - ell rows and kept when wt(e_B) <= w - t1 - t2.  (t1 = t2 = 1, 2 or 3.)
    Lee-Brickell with parameter t is the degenerate case ell = 0 (no collision); it is dominated
    by Stern here and is only kept for the probability table.
  * ``subset_mitm_iteration``: a random column subset S of size s and the exhaustive MITM on
    H[:, S] (the repo's numba kernels when available): catches a word iff supp(e) subset S.
  * ``find_low_weight(H, p, w, iters, seed, ...)``: the driver.  method = "auto" runs the
    exhaustive MITM when its list fits the budget (then iters is irrelevant and the answer is a
    proof on the column set), otherwise Stern iterations within the list budget.

Per-iteration success probability for ONE word of weight w' (w' <= w), uniform random choices:
  Stern(t1, t2, ell, z1, z2):  j = w' - t1 - t2 >= 0 positions must land in B and avoid the ell
    collision rows; the n columns are partitioned into (B: r, Z1: z1, Z2: z2, rest):
       P = [ w'! / (t1! t2! j!) ] * z1_(t1) z2_(t2) r_(j) / n_(w')  *  C(r - ell, j) / C(r, j)
    with x_(a) the falling factorial.  (Words of weight w' < t1 + t2 cannot be found.)
  Lee-Brickell(t): P = sum_{j >= w'-t} C(w', j) C(n - w', r - j) / C(n, r)  (hypergeometric).
  Subset-MITM(s): P = C(s, w') / C(n, w').
  Exhaustive MITM on all n columns: P = 1 (and absence is proved).
"Not found after N iterations" then gives, for each word that might exist, a miss probability
(1 - P)^N; with M words (the law's E_w) the chance that all are missed is ~ (1 - P)^{N M}.  This
is a MODEL statement about the randomness of the information sets, NOT a proof of absence:
nothing short of the exhaustive scan proves d > w.
"""
from __future__ import annotations
import math, time, itertools
from math import comb
import numpy as np

_HASH_MASK = (1 << 64) - 1


# ------------------------------------------------------------------------------------------------
# small F_p helpers
# ------------------------------------------------------------------------------------------------
def _inv_table(p):
    inv = np.zeros(p, dtype=np.int64)
    for x in range(1, p):
        inv[x] = pow(x, -1, p)
    return inv


def verify_word(H, p, supp, coef):
    """Exact check: H e = 0 over F_p with e supported on ``supp`` with coefficients ``coef``.
    Returns (ok, weight) where weight counts the non-zero coefficients after merging repeats."""
    H = np.asarray(H, dtype=np.int64)
    e = {}
    for c, a in zip(supp, coef):
        e[int(c)] = (e.get(int(c), 0) + int(a)) % p
    e = {c: a for c, a in e.items() if a}
    if not e:
        return False, 0
    s = np.zeros(H.shape[0], dtype=np.int64)
    for c, a in e.items():
        s = (s + a * H[:, c]) % p
    return bool(not s.any()), len(e)


def _canon_word(supp, coef, p):
    e = {}
    for c, a in zip(supp, coef):
        e[int(c)] = (e.get(int(c), 0) + int(a)) % p
    items = sorted((c, a) for c, a in e.items() if a)
    if not items:
        return [], []
    lead_inv = pow(items[0][1], -1, p)
    return [c for c, _ in items], [(a * lead_inv) % p for _, a in items]


def random_information_set(H, p, rng):
    """Random column order, greedy pivots -> (pivot columns B (len r), free columns A, Q = Binv H_A).
    Rows are the rows of the reduced matrix; row i corresponds to pivot column B[i]."""
    H = np.asarray(H, dtype=np.int64) % p
    r0, n = H.shape
    perm = rng.permutation(n)
    M = H[:, perm].copy()
    inv = _inv_table(p)
    piv_rows, piv_cols = [], []
    row = 0
    for c in range(n):
        if row == r0:
            break
        nz = np.nonzero(M[row:, c])[0]
        if len(nz) == 0:
            continue
        i = row + nz[0]
        if i != row:
            M[[row, i]] = M[[i, row]]
        M[row] = (M[row] * inv[M[row, c]]) % p
        others = np.nonzero(M[:, c])[0]
        others = others[others != row]
        if len(others):
            M[others] = (M[others] - np.outer(M[others, c], M[row])) % p
        piv_cols.append(c); row += 1
    r = row
    M = M[:r]
    free = np.array([c for c in range(n) if c not in set(piv_cols)], dtype=np.int64)
    B = perm[np.array(piv_cols, dtype=np.int64)]
    A = perm[free]
    Q = np.ascontiguousarray(M[:, free])
    return B, A, Q


# ------------------------------------------------------------------------------------------------
# pattern tables and list building
# ------------------------------------------------------------------------------------------------
def _patterns(p, t):
    """All coefficient patterns of length t with leading coefficient 1: ((p-1)^(t-1), t)."""
    if t == 1:
        return np.ones((1, 1), dtype=np.int64)
    tails = np.array(list(itertools.product(range(1, p), repeat=t - 1)), dtype=np.int64)
    return np.concatenate([np.ones((len(tails), 1), dtype=np.int64), tails], axis=1)


def _combos(z, t):
    """All t-subsets of range(z) as an (C(z,t), t) int32 array (lexicographic)."""
    if t == 1:
        return np.arange(z, dtype=np.int32)[:, None]
    if t == 2:
        I, J = np.triu_indices(z, 1)
        return np.stack([I, J], axis=1).astype(np.int32)
    out = np.fromiter(itertools.chain.from_iterable(itertools.combinations(range(z), t)),
                      dtype=np.int32, count=comb(z, t) * t)
    return out.reshape(-1, t)


def list_size(z, t, p):
    return comb(z, t) * (p - 1) ** (t - 1)


def _build_list(Q, Zcols, t, p, keyrows, R, inv, block=400_000):
    """Syndromes Q c for all weight-t normalised patterns c on columns Zcols of Q.
    Returns (keys uint64, combo index int64 (into _combos(len(Zcols), t)), pattern index int64,
    zero flag bool) where key = linear hash of the ell-part of Q c normalised by its first
    non-zero entry (projective key); zero flag marks a vanishing ell-part."""
    z = len(Zcols)
    combos = _combos(z, t)
    pats = _patterns(p, t)
    ncomb, npat = len(combos), len(pats)
    QT = np.ascontiguousarray(Q[keyrows].T)            # (z_total, ell)
    keys = np.empty(ncomb * npat, dtype=np.uint64)
    cidx = np.empty(ncomb * npat, dtype=np.int32)
    pidx = np.empty(ncomb * npat, dtype=np.int32)
    zero = np.zeros(ncomb * npat, dtype=bool)
    pos = 0
    cb = max(1, block // npat)
    for c0 in range(0, ncomb, cb):
        cm = combos[c0:c0 + cb]                          # (nc, t)
        cols = Zcols[cm]                                 # global column ids
        base = [QT[cols[:, i]] for i in range(t)]        # t arrays (nc, ell)
        for pi in range(npat):
            s = base[0].copy()
            for i in range(1, t):
                s += pats[pi, i] * base[i]
            s %= p
            nzmask = s != 0
            has = nzmask.any(axis=1)
            first = np.argmax(nzmask, axis=1)
            lead = s[np.arange(len(s)), first]
            lead[~has] = 1
            s = (s * inv[lead][:, None]) % p
            k = (s.astype(np.uint64) * R[None, :]).sum(axis=1, dtype=np.uint64)
            m = len(s)
            keys[pos:pos + m] = k
            cidx[pos:pos + m] = np.arange(c0, c0 + len(cm))
            pidx[pos:pos + m] = pi
            zero[pos:pos + m] = ~has
            pos += m
    return keys, cidx, pidx, zero, combos, pats


def _join_and_complete(Q, p, Z1, Z2, L1, L2, keyrows, restrows, wmax_rest, B, A_is_same,
                       inv, max_hits=None, block=200_000, Amap=None):
    """Join two lists on the projective key, complete each match on ``restrows`` and return the
    (e_A support/coefficients, e_B support/coefficients) candidates with wt(e_B|rest) <= wmax_rest.
    Q may be the raw H (exhaustive mode: keyrows = all rows, restrows empty, B empty)."""
    keys1, c1, p1, z1flag, combos1, pats1 = L1
    keys2, c2, p2, z2flag, combos2, pats2 = L2
    same = A_is_same and combos1.shape[1] == combos2.shape[1]
    order = np.argsort(keys1, kind="stable")
    ks = keys1[order]
    a_parts, b_parts = [], []
    chunk = 5_000_000
    for q0 in range(0, len(keys2), chunk):
        kq = keys2[q0:q0 + chunk]
        lo = np.searchsorted(ks, kq, "left"); hi = np.searchsorted(ks, kq, "right")
        cnt = hi - lo
        # same list: every entry matches itself; only keys with multiplicity >= 2 can pair distinct combos
        sel = np.nonzero(cnt >= (2 if same else 1))[0]
        if len(sel) == 0:
            continue
        cnt_s = cnt[sel]
        tot = int(cnt_s.sum())
        b_idx = np.repeat(sel + q0, cnt_s)
        offs = np.arange(tot) - np.repeat(np.cumsum(cnt_s) - cnt_s, cnt_s)
        a_idx = order[np.repeat(lo[sel], cnt_s) + offs]
        if same:
            # the two halves of a word are DISTINCT combos; see each unordered pair once
            keep = c1[a_idx] < c2[b_idx]
            a_idx = a_idx[keep]; b_idx = b_idx[keep]
        a_parts.append(a_idx); b_parts.append(b_idx)
    del order, ks
    if not a_parts:
        return []
    a_idx = np.concatenate(a_parts); b_idx = np.concatenate(b_parts); tot = len(a_idx)
    QT = Q.T
    out = []
    t1, t2 = combos1.shape[1], combos2.shape[1]
    for m0 in range(0, tot, block):
        a = a_idx[m0:m0 + block]; b = b_idx[m0:m0 + block]
        cols1 = Z1[combos1[c1[a]]]; cols2 = Z2[combos2[c2[b]]]      # (m, t)
        cf1 = pats1[p1[a]]; cf2 = pats2[p2[b]]
        s1 = np.zeros((len(a), Q.shape[0]), dtype=np.int64)
        s2 = np.zeros((len(a), Q.shape[0]), dtype=np.int64)
        for i in range(t1):
            s1 += cf1[:, i:i + 1] * QT[cols1[:, i]]
        for i in range(t2):
            s2 += cf2[:, i:i + 1] * QT[cols2[:, i]]
        s1 %= p; s2 %= p
        k1 = s1[:, keyrows]; k2 = s2[:, keyrows]
        nz1 = k1 != 0; has1 = nz1.any(axis=1)
        f1 = np.argmax(nz1, axis=1)
        # lambda with s1 + lam s2 = 0 on the key rows: lam = -k1[f] / k2[f]
        v1 = k1[np.arange(len(a)), f1]; v2 = k2[np.arange(len(a)), f1]
        lam_all = np.where(has1 & (v2 != 0), (-(v1 * inv[v2 % p])) % p, 0)
        # projective match guarantees v2 != 0 when has1 (same normalised vector); zero-key pairs
        # (both ell-parts zero) get every lambda tried
        groups = [(lam_all, np.ones(len(a), dtype=bool))] if not (~has1).any() else \
            [(lam_all, has1)] + [(np.full(len(a), lam, dtype=np.int64), ~has1) for lam in range(1, p)]
        for lam, mask in groups:
            if not mask.any():
                continue
            full = (s1 + lam[:, None] * s2) % p
            eB = (-full) % p
            if len(restrows):
                wrest = np.count_nonzero(eB[:, restrows], axis=1)
                ok = mask & (wrest <= wmax_rest)
            else:
                ok = mask & ~(full != 0).any(axis=1)
            for i in np.nonzero(ok)[0]:
                supp = list(cols1[i]) + list(cols2[i])
                if Amap is not None:
                    supp = [int(Amap[c]) for c in supp]
                coef = list(cf1[i]) + [(int(lam[i]) * int(x)) % p for x in cf2[i]]
                if len(B):
                    nzB = np.nonzero(eB[i])[0]
                    supp += [int(B[j]) for j in nzB]
                    coef += [int(eB[i, j]) for j in nzB]
                out.append((supp, coef))
                if max_hits is not None and len(out) >= max_hits:
                    return out
    return out


# ------------------------------------------------------------------------------------------------
# the three primitives
# ------------------------------------------------------------------------------------------------
def mitm_exhaustive(H, p, w, cols=None, rng=None, max_hits=1, return_all=False, block=400_000):
    """Exhaustive projective MITM for words of weight <= w (w >= 2) of ker(H) supported on
    ``cols`` (default: all columns): split w = t1 + t2 (t1 = ceil(w/2)); every word of weight
    w' in [2, w] with a t1 / t2 split is caught (words of weight < w: some split exists since
    t1 <= w' is not required... see note) -- in practice call with increasing w as the repo does.
    Returns a list of verified words [(supp, coef)] (empty = PROOF that no weight-w word with a
    t1 + t2 split exists on ``cols``; together with absence of weight < w this proves d > w)."""
    H = np.asarray(H, dtype=np.int64) % p
    r, n = H.shape
    cols = np.arange(n, dtype=np.int64) if cols is None else np.asarray(cols, dtype=np.int64)
    rng = np.random.default_rng(0) if rng is None else rng
    inv = _inv_table(p)
    R = rng.integers(1, 2 ** 63 - 1, size=r, dtype=np.int64).astype(np.uint64) | np.uint64(1)
    t1 = (w + 1) // 2; t2 = w - t1
    keyrows = np.arange(r)
    L1 = _build_list(H, cols, t1, p, keyrows, R, inv, block)
    L2 = L1 if t1 == t2 else _build_list(H, cols, t2, p, keyrows, R, inv, block)
    cands = _join_and_complete(H, p, cols, cols, L1, L2, keyrows, np.array([], dtype=np.int64), 0,
                               np.array([], dtype=np.int64), True, inv,
                               max_hits=None if return_all else max(64, 4 * max_hits))
    words = {}
    for supp, coef in cands:
        s, c = _canon_word(supp, coef, p)
        if not s or len(s) > w:
            continue
        ok, wt = verify_word(H, p, s, c)
        if ok:
            words[tuple(s)] = (s, c)
        if not return_all and len(words) >= max_hits:
            break
    return list(words.values())


def stern_iteration(H, p, w, rng, t1=2, t2=2, ell=None, z1=None, z2=None, max_hits=1, block=400_000):
    """One p-ary Stern/Dumer iteration.  Returns (verified words, info dict)."""
    H = np.asarray(H, dtype=np.int64) % p
    n = H.shape[1]
    B, A, Q = random_information_set(H, p, rng)
    r = len(B)
    nA = len(A)
    z1 = nA // 2 if z1 is None else min(z1, nA // 2)
    z2 = nA - nA // 2 if z2 is None else min(z2, nA - z1)
    permA = rng.permutation(nA)
    Z1 = permA[:z1]; Z2 = permA[z1:z1 + z2]            # positions within A (= columns of Q)
    N1, N2 = list_size(z1, t1, p), list_size(z2, t2, p)
    if ell is None:
        # about one projective key per list entry: (p^ell - 1)/(p-1) ~ 4 * max(N1, N2)
        ell = int(min(r, max(1, math.ceil(math.log(4 * max(N1, N2) * (p - 1) + 1, p)))))
    rows = rng.permutation(r)
    keyrows = np.sort(rows[:ell]); restrows = np.sort(rows[ell:])
    inv = _inv_table(p)
    R = rng.integers(1, 2 ** 63 - 1, size=ell, dtype=np.int64).astype(np.uint64) | np.uint64(1)
    L1 = _build_list(Q, Z1, t1, p, keyrows, R, inv, block)
    L2 = _build_list(Q, Z2, t2, p, keyrows, R, inv, block)
    wrest = w - t1 - t2
    cands = _join_and_complete(Q, p, Z1, Z2, L1, L2, keyrows, restrows, wrest, B, False, inv,
                               max_hits=max(64, 4 * max_hits), Amap=A)
    words = {}
    for supp, coef in cands:
        s, c = _canon_word(supp, coef, p)
        if not s or len(s) > w:
            continue
        ok, wt = verify_word(H, p, s, c)
        if ok:
            words[tuple(s)] = (s, c)
    info = {"r": r, "ell": ell, "z1": z1, "z2": z2, "N1": N1, "N2": N2, "t1": t1, "t2": t2}
    return list(words.values()), info


def subset_mitm_iteration(H, p, w, rng, s, use_repo=True):
    """Random column subset S of size s, exhaustive MITM for weight <= w on H[:, S].
    Uses the repo's numba kernels (weight_d_exists_fast / weight6_witness_3x3 /
    weight7_witness_3x4) when available and w <= 7, else the numpy MITM above."""
    H = np.asarray(H, dtype=np.int64) % p
    n = H.shape[1]
    S = np.sort(rng.choice(n, size=min(s, n), replace=False))
    Hs = np.ascontiguousarray(H[:, S])
    found = None
    if use_repo:
        try:
            from qmsd.mindist_nb import weight_d_exists_fast, weight6_witness_3x3, weight7_witness_3x4
            low = None
            for d in range(2, min(w, 5) + 1):
                if weight_d_exists_fast(Hs, p, d):
                    low = d; break
            if low is not None:
                found = mitm_exhaustive(Hs, p, low, rng=rng)
            elif w >= 6:
                wit = weight6_witness_3x3(Hs, p)
                if wit is not None:
                    found = mitm_exhaustive(Hs, p, 6, cols=np.array(wit), rng=rng)
                elif w >= 7:
                    wit = weight7_witness_3x4(Hs, p)
                    if wit is not None:
                        found = mitm_exhaustive(Hs, p, 7, cols=np.array(wit), rng=rng)
        except ImportError:
            use_repo = False
    if not use_repo:
        found = None
        for d in range(2, w + 1):
            found = mitm_exhaustive(Hs, p, d, rng=rng)
            if found:
                break
    words = []
    for supp, coef in (found or []):
        g = [int(S[c]) for c in supp]
        ok, wt = verify_word(H, p, g, coef)
        if ok:
            words.append((g, coef))
    return words, {"s": int(len(S))}


# ------------------------------------------------------------------------------------------------
# success probabilities (model statements, not proofs)
# ------------------------------------------------------------------------------------------------
def _ff(x, a):
    out = 1
    for i in range(a):
        out *= (x - i)
    return out


def p_success_stern(n, r, wprime, t1, t2, ell, z1, z2):
    j = wprime - t1 - t2
    if j < 0 or j > r - ell or t1 > z1 or t2 > z2:
        return 0.0
    multi = math.factorial(wprime) // (math.factorial(t1) * math.factorial(t2) * math.factorial(j))
    base = multi * _ff(z1, t1) * _ff(z2, t2) * _ff(r, j) / _ff(n, wprime)
    return base * comb(r - ell, j) / comb(r, j)


def p_success_lee_brickell(n, r, wprime, t):
    return sum(comb(wprime, j) * comb(n - wprime, r - j) for j in range(max(0, wprime - t), min(wprime, r) + 1)) / comb(n, r)


def p_success_subset(n, s, wprime):
    return comb(s, wprime) / comb(n, wprime) if s >= wprime else 0.0


def expected_words(n, checks, p, wprime):
    """The law: E[number of projective weight-w words] = C(n, w)(p-1)^(w-1)/p^checks."""
    return comb(n, wprime) * (p - 1) ** (wprime - 1) / p ** checks


# ------------------------------------------------------------------------------------------------
# driver
# ------------------------------------------------------------------------------------------------
def find_low_weight(H, p, w, iters=100, seed=0, method="auto", list_budget=3_000_000, t1=2, t2=2,
                    subset_size=None, time_limit=None, verbose=False, checks=None, stop_at_first=True):
    """Look for a word of ker(H) of weight <= w.  Returns a dict with keys
    found (bool), support, coefficients, weight, method, iterations, seconds, params,
    p_success (per-iteration probability of catching one word of weight w' for w' = 2..w),
    p_miss_one_word (= (1 - p_success)^iterations), model_expected_words (the law, if ``checks``
    or rank is known) and p_all_missed (~ (1 - p_success)^(iterations * E_w)).
    Every reported word is verified exactly (H e = 0 mod p, weight)."""
    H = np.ascontiguousarray(np.asarray(H, dtype=np.int64) % p)
    r0, n = H.shape
    rng = np.random.default_rng(seed)
    t_start = time.time()
    res = {"found": False, "support": None, "coefficients": None, "weight": None, "method": method,
           "iterations": 0, "seconds": 0.0, "params": {}, "p_success": {}, "n": int(n), "rows": int(r0)}
    if checks is None:
        from qmsd.structured_ad import _rank
        checks = int(_rank(H.copy(), p))
    res["checks"] = int(checks)
    res["model_expected_words"] = {str(wp): expected_words(n, checks, p, wp) for wp in range(2, w + 1)}

    def finish(words, it, extra=None):
        res["iterations"] = it; res["seconds"] = round(time.time() - t_start, 2)
        if extra:
            res["params"].update(extra)
        if words:
            s, c = min(words, key=lambda sc: len(sc[0]))
            res.update({"found": True, "support": [int(x) for x in s], "coefficients": [int(x) for x in c],
                        "weight": len(s), "all_found": [[int(x) for x in s_] for s_, _ in words]})
        for wp, ps in list(res["p_success"].items()):
            res.setdefault("p_miss_one_word", {})[wp] = (1 - ps) ** it if ps < 1 else 0.0
            E = res["model_expected_words"][wp]
            res.setdefault("p_all_missed_model", {})[wp] = math.exp(-E * (1 - (1 - ps) ** it))
        return res

    th = (w + 1) // 2
    if method == "auto":
        method = "exhaustive" if list_size(n, th, p) <= list_budget else "stern"
    res["method"] = method

    if method == "exhaustive":
        words = []
        for d in range(2, w + 1):
            words = mitm_exhaustive(H, p, d, rng=rng)
            if words:
                break
        for wp in range(2, w + 1):
            res["p_success"][str(wp)] = 1.0
        return finish(words, 1, {"list_size": list_size(n, th, p), "exhaustive": True})

    if method == "stern":
        # choose window sizes so that both lists fit the budget
        r = checks
        nA = n - r
        zmax = nA // 2
        z = zmax
        while list_size(z, max(t1, t2), p) > list_budget and z > max(t1, t2):
            z = int(z * 0.9)
        z1 = z2 = z
        words_all = []; info = None
        for it in range(1, iters + 1):
            words, info = stern_iteration(H, p, w, rng, t1=t1, t2=t2, z1=z1, z2=z2)
            if verbose:
                print(f"  stern it {it}: {len(words)} words, {time.time()-t_start:.1f}s", flush=True)
            if words:
                words_all += words
                if stop_at_first:
                    break
            if time_limit is not None and time.time() - t_start > time_limit:
                break
        for wp in range(2, w + 1):
            res["p_success"][str(wp)] = p_success_stern(n, info["r"], wp, t1, t2, info["ell"], info["z1"], info["z2"])
        return finish(words_all, it, info)

    if method == "subset":
        s = subset_size
        if s is None:
            s = n
            while list_size(s, th, p) > list_budget and s > w:
                s = int(s * 0.9)
        words_all = []
        for it in range(1, iters + 1):
            words, info = subset_mitm_iteration(H, p, w, rng, s)
            if verbose:
                print(f"  subset it {it}: {len(words)} words, {time.time()-t_start:.1f}s", flush=True)
            if words:
                words_all += words
                if stop_at_first:
                    break
            if time_limit is not None and time.time() - t_start > time_limit:
                break
        for wp in range(2, w + 1):
            res["p_success"][str(wp)] = p_success_subset(n, s, wp)
        return finish(words_all, it, {"s": s})

    if method == "lee_brickell":
        t = t1
        words_all = []
        for it in range(1, iters + 1):
            B, A, Q = random_information_set(H, p, rng)
            rr = len(B)
            words = []
            if t == 1:
                wB = np.count_nonzero(Q, axis=0)
                for a in np.nonzero(wB <= w - 1)[0]:
                    eB = (-Q[:, a]) % p; nz = np.nonzero(eB)[0]
                    words.append(([int(A[a])] + [int(B[j]) for j in nz], [1] + [int(eB[j]) for j in nz]))
            elif t == 2:
                nA = len(A)
                for a in range(nA - 1):
                    for c in range(1, p):
                        S = (Q[:, a:a + 1] + c * Q[:, a + 1:]) % p
                        wB = np.count_nonzero(S, axis=0)
                        for b in np.nonzero(wB <= w - 2)[0]:
                            eB = (-S[:, b]) % p; nz = np.nonzero(eB)[0]
                            words.append(([int(A[a]), int(A[a + 1 + b])] + [int(B[j]) for j in nz],
                                          [1, c] + [int(eB[j]) for j in nz]))
            else:
                raise ValueError("lee_brickell implemented for t = 1, 2")
            good = []
            for s_, c_ in words:
                s_, c_ = _canon_word(s_, c_, p)
                ok, wt = verify_word(H, p, s_, c_)
                if ok and wt <= w:
                    good.append((s_, c_))
            if good:
                words_all += good
                if stop_at_first:
                    break
            if time_limit is not None and time.time() - t_start > time_limit:
                break
        for wp in range(2, w + 1):
            res["p_success"][str(wp)] = p_success_lee_brickell(n, rr, wp, t)
        return finish(words_all, it, {"t": t})

    raise ValueError(method)


def cost_table(n, r, p, w, list_budget=3_000_000):
    """Expected iterations per word of weight w for each variant under the list budget (model)."""
    rows = []
    th = (w + 1) // 2
    rows.append(("exhaustive MITM (proof)", list_size(n, th, p), 1.0))
    for t in (1, 2):
        rows.append((f"Lee-Brickell t={t}", comb(n - r, t) * (p - 1) ** (t - 1), p_success_lee_brickell(n, r, w, t)))
    for (a, b) in ((2, 2), (3, 3), (1, 1), (2, 3)):
        if a + b > w:
            continue
        z = (n - r) // 2
        while list_size(z, max(a, b), p) > list_budget and z > max(a, b):
            z = int(z * 0.9)
        N = list_size(z, max(a, b), p)
        ell = int(min(r, max(1, math.ceil(math.log(4 * N * (p - 1) + 1, p)))))
        rows.append((f"Stern {a}+{b}, z={z}, ell={ell}", N, p_success_stern(n, r, w, a, b, ell, z, z)))
    s = n
    while list_size(s, th, p) > list_budget and s > w:
        s = int(s * 0.9)
    rows.append((f"subset-MITM s={s}", list_size(s, th, p), p_success_subset(n, s, w)))
    return rows


if __name__ == "__main__":
    import sys, json, argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("npz"); ap.add_argument("w", type=int)
    ap.add_argument("--method", default="auto"); ap.add_argument("--iters", type=int, default=100)
    ap.add_argument("--seed", type=int, default=0); ap.add_argument("--budget", type=int, default=3_000_000)
    ap.add_argument("--t1", type=int, default=2); ap.add_argument("--t2", type=int, default=2)
    ap.add_argument("--subset", type=int, default=None); ap.add_argument("--time", type=float, default=None)
    ap.add_argument("--out", default=None); ap.add_argument("--all", action="store_true", help="do not stop at the first word")
    ap.add_argument("--table", action="store_true", help="print the expected-iterations table and exit")
    a = ap.parse_args()
    d = np.load(a.npz, allow_pickle=True)
    H = d["H"]; p = int(d["p"]); rk = int(d["rank"])
    if a.table:
        print(f"n={H.shape[1]} r={rk} p={p} w={a.w}: E_w(law) = {expected_words(H.shape[1], rk, p, a.w):.3g}")
        for name, N, ps in cost_table(H.shape[1], rk, p, a.w, a.budget):
            print(f"  {name:36s} list {N:12.3g}  P/iter/word {ps:9.3g}  E[iters/word] {1/ps if ps else float('inf'):10.3g}")
        sys.exit(0)
    res = find_low_weight(H, p, a.w, iters=a.iters, seed=a.seed, method=a.method, list_budget=a.budget, t1=a.t1, t2=a.t2,
                          subset_size=a.subset, time_limit=a.time, verbose=True, checks=rk, stop_at_first=not a.all)
    res["instance"] = a.npz
    print(json.dumps({k: v for k, v in res.items() if k != "all_found"}, indent=1, default=float))
    if a.out:
        json.dump(res, open(a.out, "w"), indent=1, default=float)
