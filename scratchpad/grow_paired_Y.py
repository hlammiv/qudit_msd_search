"""Grow output sets Y in AG(3,7) for the paired-cap codes: randomized greedy with restarts.  A point x is added only if
Y + x keeps <= 3 points on every line, stays independent for RM(5,3) and RM(6,3), and keeps d'(Y + x) = 4 (no word of
weight <= 3 in ker of RM(5,3) shortened at Y + x).  d'' is NOT checked here (cert_d2_paired.py on lenore).
Writes the largest sets found (as probe_paired_cap.py "draws" so the cert scripts can read them).
usage: PYTHONPATH=. python scratchpad/grow_paired_Y.py RESTARTS SEED OUT.json [START_FROM_SWEEP.json T SAMPLE]"""
import sys, json, time, itertools, random, numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.triorthogonal import shorten_matrix
from qmsd.structured_ad import _rank
from qmsd.mindist_nb import weight_d_exists_fast

restarts, seed, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]; p = 7
P3 = np.asarray(points(3, p), dtype=np.int64); idx3 = {tuple(int(v) for v in x): i for i, x in enumerate(P3)}
R5 = np.asarray(rm_generator(5, 3, p), dtype=np.int64) % p; R6 = np.asarray(rm_generator(6, 3, p), dtype=np.int64) % p
_dirs = sorted({min(tuple(int(v) for v in (c * np.array(d)) % p) for c in range(1, p)) for d in itertools.product(range(p), repeat=3) if any(d)})
LINES = sorted({tuple(sorted(idx3[tuple(int(v) for v in (P3[a] + c * np.array(d)) % p)] for c in range(p))) for a in range(p ** 3) for d in _dirs})
THROUGH = [[] for _ in range(p ** 3)]
for li, L in enumerate(LINES):
    for x in L: THROUGH[x].append(li)
start_Y = None
if len(sys.argv) > 4:
    start_Y = [d for d in json.load(open(sys.argv[4]))["draws"] if d["t"] == int(sys.argv[5]) and d["sample"] == int(sys.argv[6])][0]["Y"]

def dprime_ge4(Y):
    H = np.asarray(shorten_matrix(R5, Y, p), dtype=np.int64) % p
    return not weight_d_exists_fast(H, p, 2) and not weight_d_exists_fast(H, p, 3)

SWAPS = int(__import__('os').environ.get('SWAPS', '0'))
rng = random.Random(seed); best = []; T0 = time.time(); draws = []
for rs in range(restarts):
    if start_Y is not None and rs == 0:
        Y = list(start_Y)
    else:                                          # fresh random no-4-collinear seed of size 36 (d' = 4 is automatic-ish there)
        Y = []
        order = list(range(p ** 3)); rng.shuffle(order); cnt = np.zeros(len(LINES), dtype=np.int64)
        for x in order:
            if (cnt[THROUGH[x]] < 3).all():
                Y.append(x); cnt[THROUGH[x]] += 1
                if len(Y) == 36: break
        if not dprime_ge4(Y) or _rank(R5[:, Y], p) < len(Y): continue
    cnt = np.zeros(len(LINES), dtype=np.int64)
    for x in Y: cnt[THROUGH[x]] += 1
    while True:
        cands = [x for x in range(p ** 3) if x not in set(Y) and (cnt[THROUGH[x]] < 3).all()]; rng.shuffle(cands); added = False
        for x in cands:
            Z = Y + [x]
            if _rank(R5[:, Z], p) < len(Z) or _rank(R6[:, Z], p) < len(Z): continue
            if dprime_ge4(Z):
                Y = Z; cnt[THROUGH[x]] += 1; added = True; break
        if not added:                              # swap move: drop one point, try to add two (up to SWAPS tries)
            improved = False
            for _ in range(SWAPS):
                drop = rng.choice(Y); Yd = [y for y in Y if y != drop]
                cd = np.zeros(len(LINES), dtype=np.int64)
                for y in Yd: cd[THROUGH[y]] += 1
                c2 = [x for x in range(p ** 3) if x not in set(Y) and (cd[THROUGH[x]] < 3).all()]; rng.shuffle(c2)
                got = []
                for x in c2:
                    Z = Yd + got + [x]
                    cz = cd.copy()
                    for y in got: cz[THROUGH[y]] += 1
                    if not (cz[THROUGH[x]] < 3).all(): continue
                    if _rank(R5[:, Z], p) < len(Z) or _rank(R6[:, Z], p) < len(Z) or not dprime_ge4(Z): continue
                    got.append(x)
                    if len(got) == 2: break
                if len(got) == 2:
                    Y = Yd + got; cnt = np.zeros(len(LINES), dtype=np.int64)
                    for y in Y: cnt[THROUGH[y]] += 1
                    improved = True; break
            if not improved: break
    n, k = 686 - 2 * len(Y), 2 * len(Y)
    print(f"[{time.time()-T0:5.0f}s] restart {rs}: maximal |Y| = {len(Y)} -> [[{n},{k},8 if d''>=8]] gamma {np.log(n/k)/np.log(8):.3f}", flush=True)
    draws.append({"t": len(Y), "sample": rs, "Y": sorted(Y), "d1": 4, "seed": seed})
    json.dump({"construction": "greedy grown Y (grow_paired_Y.py), d' = 4, d'' unchecked", "draws": draws}, open(out, "w"), indent=1)
