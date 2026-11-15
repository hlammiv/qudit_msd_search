"""Paired-cap outputs on the single-species support g = x1(x1-1)...(x1-4) (p = 7, m = 4, r = 6; X = {x1 = 5} u {x1 = 6},
charges +-1, one Clifford species).  RM(6,4)|X = {(a, a+b): a in RM(6,3), b in RM(5,3)}; with outputs S = Y x {5,6}
the distance is exactly d = min(2 d'(Y), d''(Y)), where d' (d'') is the minimum number of points off Y over nonzero
words of RM(12,3) = RM(5,3)^perp (RM(11,3) = RM(6,3)^perp), i.e. the distance of the m = 3 code shortened at Y.
Y = random cap of AG(3,7) (no 3 collinear) of size t, independent for RM(6,3) and RM(5,3), so k = 2t, n = 686 - 2t.
This script: d' exactly (weights <= 5) and d'' up to weight W2 (default 5); full-code cross-check to weight 4.
usage: PYTHONPATH=. python scratchpad/probe_paired_cap.py T1,T2,.. DRAWS SEED OUT.json [W2 [MAXLINE=3]]"""
import sys, random, time, json, itertools, numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.weighted import weighted_support, product_of_linear
from qmsd.structured_ad import _rank
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import weight_d_exists_fast, weight5_exists_parallel

ts = [int(v) for v in sys.argv[1].split(",")]; draws, seed, out_path = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
W2 = int(sys.argv[5]) if len(sys.argv) > 5 else 5
MAXLINE = int(sys.argv[6]) if len(sys.argv) > 6 else 3
p = 7
P3 = np.asarray(points(3, p), dtype=np.int64); idx3 = {tuple(int(v) for v in x): i for i, x in enumerate(P3)}
R5 = np.asarray(rm_generator(5, 3, p), dtype=np.int64) % p; R6 = np.asarray(rm_generator(6, 3, p), dtype=np.int64) % p
g = product_of_linear(p, 4, [((1, 0, 0, 0), c) for c in range(5)])
GX, X, mu = weighted_support(p, 4, g, 6)
P4 = np.asarray(points(4, p), dtype=np.int64)[X]; idxX = {tuple(int(v) for v in x): i for i, x in enumerate(P4)}

# all lines of AG(3,7) as point-index arrays, and the lines through each point
_dirs = sorted({min(tuple(int(v) for v in (c * np.array(d)) % p) for c in range(1, p))
                for d in itertools.product(range(p), repeat=3) if any(d)})
LINES = sorted({tuple(sorted(idx3[tuple(int(v) for v in (P3[a] + c * np.array(d)) % p)] for c in range(p)))
                for a in range(p ** 3) for d in _dirs})
THROUGH = [[] for _ in range(p ** 3)]
for li, L in enumerate(LINES):
    for x in L: THROUGH[x].append(li)
assert len(LINES) == 2793

def random_cap(t, rng, maxline=2, tries=200):
    """Random greedy set with at most ``maxline`` points on any line (maxline = 2: a cap), independent for RM(5,3)
    and RM(6,3)."""
    for _ in range(tries):
        Y = []; cnt = np.zeros(len(LINES), dtype=np.int64); order = list(range(p ** 3)); rng.shuffle(order)
        for x in order:
            if (cnt[THROUGH[x]] < maxline).all():
                Y.append(x); cnt[THROUGH[x]] += 1
                if len(Y) == t: break
        if len(Y) == t and _rank(R5[:, Y], p) == t and _rank(R6[:, Y], p) == t: return sorted(Y)
    raise RuntimeError(f"no independent set of size {t} (maxline {maxline}) in {tries} tries")

def min_off(R, Y, wmax):
    """min weight of ker(shorten(R, Y)) up to wmax (None if > wmax)."""
    H = np.asarray(shorten_matrix(R, Y, p), dtype=np.int64) % p
    for w in range(2, min(wmax, 4) + 1):
        if weight_d_exists_fast(H, p, w): return w
    if wmax >= 5 and weight5_exists_parallel(H, p): return 5
    return None

rng = random.Random(seed); out = []
for t in ts:
    for s in range(draws):
        T0 = time.time(); Y = random_cap(t, rng, MAXLINE)
        d1 = min_off(R5, Y, 5); d2 = min_off(R6, Y, W2)
        S = sorted(idxX[(c,) + tuple(int(v) for v in P3[y])] for y in Y for c in (5, 6))
        assert _rank(GX[:, S].T.copy(), p) == 2 * t
        G0 = np.asarray(shorten_matrix(GX, S, p), dtype=np.int64) % p
        dfull = next((w for w in (2, 3, 4) if weight_d_exists_fast(G0, p, w)), None)
        lo = min(2 * (d1 if d1 else 6), d2 if d2 else W2 + 1)
        n, k = 686 - 2 * t, 2 * t
        print(f"t={t} s={s}: [[{n},{k}]] d'={d1 if d1 else '>=6'} d''={d2 if d2 else '>' + str(W2)} "
              f"-> d=min(2d',d'') {'=' if d2 else '>='}{lo}; gamma{'=' if d2 else '<='}{np.log(n/k)/np.log(lo):.3f}; "
              f"full-code weight<=4 word: {dfull} [{time.time()-T0:.0f}s]", flush=True)
        assert dfull is None or dfull >= lo or lo > 4
        out.append({"maxline": MAXLINE, "t": t, "sample": s, "n": n, "k": k, "Y": Y, "d1": d1, "d2": d2, "d2_checked_to": W2, "d_lower": lo, "full_w4": dfull})
        json.dump({"construction": "S = Y x {5,6}, Y cap in AG(3,7)", "seed": seed, "draws": out}, open(out_path, "w"), indent=1)
