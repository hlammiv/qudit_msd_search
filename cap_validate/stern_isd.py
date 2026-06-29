"""Optimised multithreaded q-ary (F_3) Stern / ISD low-weight-codeword finder.

METHOD 1 of the cap-code validator.  Given a parity-check matrix ``H`` (= the
X-stabiliser ``G0``, an ``r x n`` matrix over F_3 with ``r`` small, ``n`` large) we
search for low-weight nonzero codewords ``c`` of the dual code ``C = G0^perp =
null(H)`` -- equivalently the minimum number of F_3-linearly-dependent columns of
``H``.  For the cap code ``H`` is ``55 x 1968`` and the validation question is
whether ``d(C) >= 10``.

Algorithm (q-ary Stern, the "find a low-weight codeword" / homogeneous variant):

  1. Random column permutation; F_3 Gaussian elimination of ``H`` to RREF.  This
     selects ``r`` *pivot* columns (forming I_r) and ``k = n - r`` *non-pivot*
     columns.  After RREF the non-pivot column ``c`` is its own coordinate vector
     ``CV[c] in F_3^r`` w.r.t. the pivot basis.  A codeword is obtained by choosing
     values on the non-pivot (information) coordinates; the pivot coordinates are
     then forced to ``-sum_c v_c CV[c]``.  Codeword weight = (#nonzero info coords)
     + (#nonzero pivot coords).

  2. Split the ``k`` non-pivot coords into two halves A, B.  Pick ``l`` of the ``r``
     pivot rows as a collision "window".  Enumerate all size-``p`` column subsets of
     A (each chosen column with an F_3 coeff in {1,2}); key each by its window
     vector.  Enumerate size-``p`` subsets of B keyed by the *negated* window vector.
     A collision (matching window key) means the combined window coordinates cancel
     -- a codeword candidate supported on ``2p`` info coords.  Compute its full pivot
     weight (with early-exit pruning) and record total = ``2p + pivot_weight``.

  3. Repeat over many random permutations / windows (deep search), keeping the global
     lowest weight and a witness codeword.

One-sided guarantee: any codeword of weight ``w < 10`` that is *found* PROVES
``d < 10`` (and is returned as an explicit witness); finding nothing ``< 10`` after a
deep search is strong probabilistic evidence (not proof) that ``d >= 10``.

The heavy per-iteration work (RREF, sublist construction, birthday collision, weight
evaluation) is a numba ``nogil`` njit kernel, so a ThreadPoolExecutor parallelises
iterations across cores with the GIL released.
"""
from __future__ import annotations

import math
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

import numpy as np
from numba import njit

# F_3 multiplicative inverse table: inv[1]=1, inv[2]=2 (2*2=4=1 mod 3).
_INV3 = np.array([0, 1, 2], dtype=np.int8)


# ---------------------------------------------------------------------------
# Numba kernels
# ---------------------------------------------------------------------------
@njit(cache=True, nogil=True)
def _rref_f3(Hw, r, n):
    """In-place F_3 RREF of ``Hw`` (r x n, int8, already column-permuted).

    Returns (pivot_pos, nonpivot_pos): pivot_pos[a] = permuted-column position whose
    column became the a-th unit vector; nonpivot_pos lists the remaining columns in
    increasing position order.  Assumes ``H`` has rank r (full row rank); pivots are
    taken in increasing column position so a random permutation gives a random info
    set.
    """
    pivot_pos = np.empty(r, np.int64)
    is_pivot = np.zeros(n, np.uint8)
    prow = 0
    for col in range(n):
        if prow >= r:
            break
        # find a pivot in this column at or below prow
        piv = -1
        for row in range(prow, r):
            if Hw[row, col] != 0:
                piv = row
                break
        if piv == -1:
            continue
        # swap piv <-> prow
        if piv != prow:
            for cc in range(n):
                tmp = Hw[prow, cc]
                Hw[prow, cc] = Hw[piv, cc]
                Hw[piv, cc] = tmp
        # normalise pivot row so leading entry == 1
        lead = Hw[prow, col]
        if lead == 2:  # inv(2)=2
            for cc in range(n):
                v = Hw[prow, cc]
                if v != 0:
                    Hw[prow, cc] = (v * 2) % 3
        # eliminate this column from all other rows
        for row in range(r):
            if row == prow:
                continue
            f = Hw[row, col]
            if f != 0:
                # row <- row - f * prow  (mod 3)
                for cc in range(n):
                    pv = Hw[prow, cc]
                    if pv != 0:
                        Hw[row, cc] = (Hw[row, cc] - f * pv) % 3
        pivot_pos[prow] = col
        is_pivot[col] = 1
        prow += 1

    nonpivot_pos = np.empty(n - r, np.int64)
    idx = 0
    for col in range(n):
        if is_pivot[col] == 0:
            nonpivot_pos[idx] = col
            idx += 1
    return pivot_pos, nonpivot_pos


@njit(cache=True, nogil=True)
def _comb(nh, p):
    """Binomial C(nh, p) for small p (numba has no math.comb)."""
    if p < 0 or p > nh:
        return 0
    num = 1
    den = 1
    for i in range(p):
        num *= nh - i
        den *= i + 1
    return num // den


@njit(cache=True, nogil=True)
def _comb_first(p):
    a = np.empty(p, np.int64)
    for i in range(p):
        a[i] = i
    return a


@njit(cache=True, nogil=True)
def _comb_next(a, p, nh):
    """Advance combination odometer ``a`` (p indices into 0..nh-1). False when done."""
    i = p - 1
    while i >= 0:
        if a[i] != i + nh - p:
            a[i] += 1
            for j in range(i + 1, p):
                a[j] = a[j - 1] + 1
            return True
        i -= 1
    return False


@njit(cache=True, nogil=True)
def _build_sublist(CVw, half, p, neg, pow3, out_key, out_cols, out_coef):
    """Fill a Stern sublist: every (size-p column subset of ``half``, coeff pattern).

    CVw: (k, l) window slice of CV (int8).  half: int64 indices into CV.
    Stores window key (int64), the p CV-indices, and the p coeffs (in {1,2}).
    If ``neg`` the key encodes the *negated* window vector (mod 3).
    Returns the number of entries written.
    """
    nh = half.shape[0]
    l = CVw.shape[1]
    cnt = 0
    if p > nh:
        return 0
    a = _comb_first(p)
    npat = 1 << p
    while True:
        for pat in range(npat):
            key = 0
            for t in range(l):
                s = 0
                for j in range(p):
                    cf = 1 + ((pat >> j) & 1)  # 1 or 2
                    s += cf * CVw[half[a[j]], t]
                s %= 3
                if neg and s != 0:
                    s = 3 - s
                key += s * pow3[t]
            out_key[cnt] = key
            base = cnt * p
            for j in range(p):
                out_cols[base + j] = half[a[j]]
                out_coef[base + j] = 1 + ((pat >> j) & 1)
            cnt += 1
        if not _comb_next(a, p, nh):
            break
    return cnt


@njit(cache=True, nogil=True)
def _eval_candidate(CV, r, cols, coef, twop, pivot_budget):
    """Pivot-coordinate weight of a candidate (early-exit at >pivot_budget).

    Returns (pivot_weight, acc) where acc[a] = (sum coef*CV[col])[a] mod 3.
    pivot_weight = -1 if it exceeds pivot_budget (pruned).
    """
    acc = np.zeros(r, np.int8)
    for j in range(twop):
        c = cols[j]
        f = coef[j]
        for a in range(r):
            v = CV[c, a]
            if v != 0:
                acc[a] = (acc[a] + f * v) % 3
    w = 0
    for a in range(r):
        if acc[a] != 0:
            w += 1
            if w > pivot_budget:
                return -1, acc
    return w, acc


@njit(cache=True, nogil=True)
def _stern_once(H, r, n, p, l, seed, weight_budget, info_width, witness):
    """One Stern iteration.  Returns the best total weight found (>=0) or -1 if none
    within ``weight_budget``.  On improvement, ``witness`` (len n, int8) is filled
    with the codeword in ORIGINAL column order.
    """
    np.random.seed(seed)
    k = n - r

    # random column permutation
    perm = np.arange(n)
    for i in range(n - 1, 0, -1):
        j = np.random.randint(0, i + 1)
        tmp = perm[i]
        perm[i] = perm[j]
        perm[j] = tmp

    # permuted working copy
    Hw = np.empty((r, n), np.int8)
    for i in range(r):
        for jj in range(n):
            Hw[i, jj] = H[i, perm[jj]] % 3

    pivot_pos, nonpivot_pos = _rref_f3(Hw, r, n)

    # CV[c, a] = coordinate vector of non-pivot column c (index 0..k-1)
    CV = np.empty((k, r), np.int8)
    for c in range(k):
        pos = nonpivot_pos[c]
        for a in range(r):
            CV[c, a] = Hw[a, pos]

    # random window of l pivot rows
    rowperm = np.arange(r)
    for i in range(r - 1, 0, -1):
        j = np.random.randint(0, i + 1)
        tmp = rowperm[i]
        rowperm[i] = rowperm[j]
        rowperm[j] = tmp
    window = rowperm[:l].copy()
    pow3 = np.empty(l, np.int64)
    pw = 1
    for t in range(l):
        pow3[t] = pw
        pw *= 3
    CVw = np.empty((k, l), np.int8)
    for c in range(k):
        for t in range(l):
            CVw[c, t] = CV[c, window[t]]

    # choose the info coordinates to enumerate over (all k, or a random width-W
    # subset to keep the p=3 sublists tractable -- coverage restored by randomising
    # the subset across iterations), then split into two halves.
    sel = np.arange(k)
    if info_width > 0 and info_width < k:
        W = info_width
        for i in range(W):
            j = np.random.randint(i, k)
            tmp = sel[i]
            sel[i] = sel[j]
            sel[j] = tmp
    else:
        W = k
    k1 = W // 2
    halfA = sel[:k1].copy()
    halfB = sel[k1:W].copy()

    capA = _comb(k1, p) * (1 << p)
    capB = _comb(W - k1, p) * (1 << p)

    keyA = np.empty(capA, np.int64)
    colsA = np.empty(capA * p, np.int64)
    coefA = np.empty(capA * p, np.int8)
    keyB = np.empty(capB, np.int64)
    colsB = np.empty(capB * p, np.int64)
    coefB = np.empty(capB * p, np.int8)

    nA = _build_sublist(CVw, halfA, p, False, pow3, keyA, colsA, coefA)
    nB = _build_sublist(CVw, halfB, p, True, pow3, keyB, colsB, coefB)

    ordA = np.argsort(keyA[:nA], kind="mergesort")
    ordB = np.argsort(keyB[:nB], kind="mergesort")

    twop = 2 * p
    pivot_budget = weight_budget - twop
    best = -1

    ia = 0
    ib = 0
    cand_cols = np.empty(twop, np.int64)
    cand_coef = np.empty(twop, np.int8)
    while ia < nA and ib < nB:
        ka = keyA[ordA[ia]]
        kb = keyB[ordB[ib]]
        if ka < kb:
            ia += 1
        elif kb < ka:
            ib += 1
        else:
            # equal-key groups [ia, ja) x [ib, jb)
            ja = ia
            while ja < nA and keyA[ordA[ja]] == ka:
                ja += 1
            jb = ib
            while jb < nB and keyB[ordB[jb]] == kb:
                jb += 1
            for ua in range(ia, ja):
                if pivot_budget < 0:
                    break  # best == 2p already, unbeatable with this p
                ea = ordA[ua]
                baseA = ea * p
                for j in range(p):
                    cand_cols[j] = colsA[baseA + j]
                    cand_coef[j] = coefA[baseA + j]
                for ub in range(ib, jb):
                    eb = ordB[ub]
                    baseB = eb * p
                    for j in range(p):
                        cand_cols[p + j] = colsB[baseB + j]
                        cand_coef[p + j] = coefB[baseB + j]
                    pw_, acc = _eval_candidate(
                        CV, r, cand_cols, cand_coef, twop, pivot_budget
                    )
                    if pw_ < 0:
                        continue
                    total = twop + pw_
                    if best == -1 or total < best:
                        best = total
                        # reconstruct witness in original column order
                        for jj in range(n):
                            witness[jj] = 0
                        for j in range(twop):
                            c = cand_cols[j]
                            witness[perm[nonpivot_pos[c]]] = cand_coef[j]
                        for a in range(r):
                            val = acc[a]
                            if val != 0:
                                witness[perm[pivot_pos[a]]] = (3 - val) % 3
                        # tighten budget so later candidates must beat it
                        pivot_budget = best - twop - 1
                        if pivot_budget < 0:
                            break
            ia = ja
            ib = jb
    return best


# ---------------------------------------------------------------------------
# Python driver
# ---------------------------------------------------------------------------
@dataclass
class SternResult:
    best_weight: int            # lowest codeword weight found (-1 if none <= budget)
    witness: np.ndarray | None  # codeword (len n, F_3) achieving best_weight
    iterations: int
    elapsed: float
    p: int
    l: int
    n: int
    r: int
    weight_budget: int
    history: list = field(default_factory=list)  # (iter, weight) improvements

    def verify(self, H: np.ndarray) -> bool:
        """Check the witness really is a nonzero codeword of weight best_weight."""
        if self.witness is None:
            return self.best_weight == -1
        c = (self.witness.astype(np.int64)) % 3
        if not np.any(c):
            return False
        syndrome = (np.asarray(H, dtype=np.int64) @ c) % 3
        if np.count_nonzero(syndrome) != 0:
            return False
        return int(np.count_nonzero(c)) == self.best_weight


def stern_search(
    H,
    p: int = 2,
    l: int = 13,
    iterations: int = 2000,
    weight_budget: int = 14,
    threads: int = 1,
    seed: int = 0,
    target: int | None = 10,
    info_width: int = 0,
    max_list_bytes: float = 6e9,
    progress_every: float = 0.0,
    verbose: bool = False,
) -> SternResult:
    """Run the multithreaded Stern low-weight-codeword search on parity check ``H``.

    H            : r x n matrix over F_3 (the X-stabiliser G0); rank must be r.
    p            : Stern sublist subset size per half (2 is the workhorse; total info
                   support = 2p).  p in {1,2,3}; p>=3 explodes the sublists unless
                   ``info_width`` restricts the enumeration set.
    l            : collision-window size (number of pivot rows).
    iterations   : number of random-permutation iterations (split across threads).
    weight_budget: only codewords of total weight <= this are recorded.
    threads      : worker threads (numba kernel releases the GIL).
    seed         : base RNG seed; iteration i uses seed + i.
    target       : stop early if a codeword of weight < target is found (None = never).
    info_width   : if >0, enumerate from only this many (randomly chosen each
                   iteration) of the k information coords -- needed to keep p=3
                   tractable.  0 = use all k.
    max_list_bytes: per-sublist memory guard; raises if exceeded.
    """
    H = np.asarray(H, dtype=np.int8) % 3
    r, n = H.shape
    rank = int(np.linalg.matrix_rank(H.astype(np.int64)))  # sanity; small r
    if rank != r:
        raise ValueError(f"H must have full row rank r={r}; got rank {rank}")
    if p < 1:
        raise ValueError("p must be >= 1")

    # memory guard: estimate the larger sublist size and its byte cost
    k = n - r
    W = info_width if (0 < info_width < k) else k
    k1 = W // 2
    cap = max(math.comb(k1, p), math.comb(W - k1, p)) * (1 << p)
    bytes_per_entry = 8 + (8 + 1) * p + 16  # key + (col,coef)*p + argsort overhead
    est = cap * bytes_per_entry
    if est > max_list_bytes:
        raise ValueError(
            f"Stern sublist too large for p={p}, info_width={W} (~{est/1e9:.1f} GB "
            f"> {max_list_bytes/1e9:.1f} GB). Lower p or set a smaller info_width."
        )

    best_lock_weight = [weight_budget + 1]
    best_witness = [None]
    history = []
    done = [0]
    start = time.time()
    stop_flag = [False]
    last_report = [start]

    def worker(it_indices):
        local_witness = np.zeros(n, np.int8)
        local_best = weight_budget + 1
        local_wit = None
        for it in it_indices:
            if stop_flag[0]:
                break
            w = _stern_once(
                H, r, n, p, l, seed + it, weight_budget, info_width, local_witness
            )
            done[0] += 1
            if w != -1 and w < local_best:
                local_best = w
                local_wit = local_witness.copy()
                # publish
                if w < best_lock_weight[0]:
                    best_lock_weight[0] = w
                    best_witness[0] = local_wit.copy()
                    history.append((done[0], w))
                    if verbose:
                        print(f"  [iter ~{done[0]}] new best weight = {w}", flush=True)
                    if target is not None and w < target:
                        stop_flag[0] = True
            if verbose and progress_every > 0:
                now = time.time()
                if now - last_report[0] >= progress_every:
                    last_report[0] = now
                    print(
                        f"  ... {done[0]}/{iterations} iters, "
                        f"best={best_lock_weight[0] if best_lock_weight[0] <= weight_budget else None}, "
                        f"{now - start:.0f}s",
                        flush=True,
                    )
        return local_best

    # split iteration indices round-robin across threads
    all_iters = list(range(iterations))
    chunks = [all_iters[t::threads] for t in range(threads)]
    if threads == 1:
        worker(all_iters)
    else:
        with ThreadPoolExecutor(max_workers=threads) as ex:
            list(ex.map(worker, chunks))

    elapsed = time.time() - start
    bw = best_lock_weight[0]
    if bw > weight_budget:
        bw = -1
        wit = None
    else:
        wit = best_witness[0]
    return SternResult(
        best_weight=bw,
        witness=wit,
        iterations=done[0],
        elapsed=elapsed,
        p=p,
        l=l,
        n=n,
        r=r,
        weight_budget=weight_budget,
        history=history,
    )


# ---------------------------------------------------------------------------
# Cap-code loader
# ---------------------------------------------------------------------------
def stern_search_multi(
    H,
    p_values=(1, 2, 3),
    l: int = 13,
    iterations: int = 2000,
    weight_budget: int = 14,
    threads: int = 1,
    seed: int = 0,
    target: int | None = 10,
    info_width=None,
    verbose: bool = False,
) -> SternResult:
    """Run :func:`stern_search` for several sublist sizes ``p`` and merge the results.

    A Stern iteration with parameter ``p`` only finds codewords whose information-set
    support is exactly ``2p`` (``p`` in each half), so to cover all low weights one must
    sweep ``p``: ``p=1`` catches weight-w codewords with 2 info coords (e.g. weight 3),
    ``p=2`` with 4, ``p=3`` with 6, etc.  ``iterations`` is split evenly across the
    ``p`` values.  Returns the single global best (lowest weight + witness).
    """
    H = np.asarray(H, dtype=np.int8) % 3
    per = max(1, iterations // len(p_values))
    best = None
    total_iters = 0
    # info_width per p: dict {p: width}, a scalar, or None (auto: full for p<=2)
    for p in p_values:
        lp = min(l, H.shape[0])
        if isinstance(info_width, dict):
            iw = info_width.get(p, 0)
        elif info_width is None:
            iw = 0  # full enumeration (fine for p in {1,2})
        else:
            iw = int(info_width)
        res = stern_search(
            H, p=p, l=lp, iterations=per, weight_budget=weight_budget,
            threads=threads, seed=seed + 10_000_000 * p, target=target,
            info_width=iw, verbose=verbose,
        )
        total_iters += res.iterations
        if verbose:
            print(f"  [p={p}] best={res.best_weight} ({res.iterations} iters)", flush=True)
        better = res.best_weight != -1 and (
            best is None or best.best_weight == -1 or res.best_weight < best.best_weight
        )
        if best is None or better:
            best = res
        if target is not None and res.best_weight != -1 and res.best_weight < target:
            best.iterations = total_iters
            return best
    if best is not None:
        best.iterations = total_iters
    return best


def trivial_low_weight(H) -> dict:
    """Report trivial low-weight codewords of ``H^perp`` directly from the columns.

    - a zero column j  => unit vector e_j is a weight-1 codeword of H^perp (a single
      F_3-linearly-dependent column);
    - two F_3-proportional columns i,j => a weight-2 codeword (e_i - c e_j).

    These are exact, deterministic facts about H requiring no search; the Stern engine
    will also find them, but this gives an O(n) up-front summary.
    """
    H = np.asarray(H, dtype=np.int64) % 3
    r, n = H.shape
    nz = np.count_nonzero(H, axis=0)
    zero_cols = np.nonzero(nz == 0)[0]
    # normalise each nonzero column (scale so first nonzero entry == 1) and bucket
    seen = {}
    prop_pair = None
    for j in range(n):
        if nz[j] == 0:
            continue
        col = H[:, j]
        lead = col[np.nonzero(col)[0][0]]
        key = tuple((col * lead) % 3)  # inv(1)=1, inv(2)=2 in F_3
        if key in seen:
            prop_pair = (seen[key], j)
        else:
            seen[key] = j
    min_w = 1 if len(zero_cols) > 0 else (2 if prop_pair is not None else None)
    return {
        "n": n,
        "r": r,
        "num_zero_columns": int(len(zero_cols)),
        "zero_columns": zero_cols.tolist(),
        "has_proportional_pair": prop_pair is not None,
        "proportional_pair": prop_pair,
        "trivial_min_weight": min_w,
    }


def build_cap_g0(json_path: str = "cap_qutrit_code.json"):
    """Build the cap qutrit code's X-stabiliser G0 (55 x 1968 over F_3)."""
    import json

    from qmsd import reedmuller, triorthogonal

    S = json.load(open(json_path))["puncture_columns_1indexed"]
    G = reedmuller.rm_generator(4, 7, 3)
    built = triorthogonal.build_triorthogonal_code(3, 7, 4, S, G=G)
    G0 = np.asarray(built["X_stab"], dtype=np.int8) % 3
    return G0


if __name__ == "__main__":
    import argparse
    import os
    import sys

    # allow `python cap_validate/stern_isd.py ...` from the project root to import qmsd
    sys.path.insert(0, os.getcwd())

    ap = argparse.ArgumentParser(description="Stern/ISD low-weight-codeword finder")
    ap.add_argument("--json", default="cap_qutrit_code.json")
    ap.add_argument("-p", type=int, default=2)
    ap.add_argument("-l", type=int, default=13)
    ap.add_argument("--iters", type=int, default=2000)
    ap.add_argument("--budget", type=int, default=14)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--target", type=int, default=10)
    args = ap.parse_args()

    G0 = build_cap_g0(args.json)
    print(f"G0 shape {G0.shape}")
    triv = trivial_low_weight(G0)
    print(
        f"trivial structure: {triv['num_zero_columns']} zero columns, "
        f"proportional pair={triv['proportional_pair']}, "
        f"trivial_min_weight={triv['trivial_min_weight']}"
    )
    res = stern_search(
        G0,
        p=args.p,
        l=args.l,
        iterations=args.iters,
        weight_budget=args.budget,
        threads=args.threads,
        seed=args.seed,
        target=args.target,
        verbose=True,
        progress_every=10.0,
    )
    print(f"best_weight={res.best_weight}  iters={res.iterations}  {res.elapsed:.1f}s")
    print(f"witness valid: {res.verify(G0)}")
    if res.witness is not None:
        print("support:", np.nonzero(res.witness)[0].tolist())
