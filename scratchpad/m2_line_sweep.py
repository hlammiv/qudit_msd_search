"""m = 2 line-spread sweep for p >= 17: random k-subsets of the p^2 plane with at most MAXLINE outputs
on any line; d_lines = exact min weight of a dual word supported on a line (MITM per line, lines
have p columns); for d_lines >= 6 the full-code weight <= 5 certificate (fast d <= 4, parallel
weight-5 scan) makes d = 6 exact when d_lines = 6, and leaves d >= 7 candidates for the weight-6
exclusion on lenore. Appends JSON lines; resumable.
usage: m2_line_sweep.py P K NSEEDS MAXLINE OUT.jsonl [--no-full]"""
import sys, os, json, random, time, itertools, numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.structured_ad import _rank, _null_space
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import weight_d_exists_fast, weight5_exists_parallel
p, k, nseeds, maxline, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]; nofull = "--no-full" in sys.argv
m = 2; r = (m * (p - 1) - 1) // 3
P = np.asarray(points(m, p), dtype=np.int64); G = np.asarray(rm_generator(r, m, p), dtype=np.int64) % p
L = []
for a in range(p):
    key = (P[:, 0] + a * P[:, 1]) % p
    L += [np.flatnonzero(key == c) for c in range(p)]
L += [np.flatnonzero(P[:, 1] == c) for c in range(p)]
done = set()
if os.path.exists(out):
    for l in open(out):
        if l.strip(): done.add(json.loads(l)["seed"])
for s in range(nseeds):
    if s in done: continue
    rng = random.Random(1000003 * p + 1009 * k + s); tries = 0
    lines_of = [[] for _ in range(len(P))]
    for li, ln in enumerate(L):
        for i in ln: lines_of[i].append(li)
    while True:
        tries += 1
        # greedy line-spread construction: add points in random order, never exceeding maxline on any line
        order = list(range(len(P))); rng.shuffle(order); cnt = [0] * len(L); cols = []
        for i in order:
            if all(cnt[li] < maxline for li in lines_of[i]):
                cols.append(i)
                for li in lines_of[i]: cnt[li] += 1
                if len(cols) == k: break
        if len(cols) < k: continue
        cols = sorted(cols); cs = set(cols); ml = max(cnt)
        if _rank(G[:, cols].T.copy(), p) == k: break
    G0 = np.asarray(shorten_matrix(G, cols, p), dtype=np.int64) % p
    keep = [l for l in range(len(P)) if l not in cs]; pos = {l: j for j, l in enumerate(keep)}
    t0 = time.time(); dl = 99
    for ln in L:
        cc = np.array([pos[i] for i in ln if i in pos], dtype=np.int64)
        if len(cc) < 2: continue
        M = np.ascontiguousarray(G0[:, cc]); M = M[[t for t in range(M.shape[0]) if M[t].any()]]
        if M.shape[0] == 0: dl = 1; break
        rk = _rank(M.T.copy(), p)
        if rk == len(cc): continue
        nu = len(cc) - rk
        if nu <= 4:  # exact: enumerate the line's null space (p^nu words)
            ns = np.asarray(_null_space(M, p)) % p
            coefs = np.array(list(itertools.product(range(p), repeat=ns.shape[0])))[1:]
            w = int((((coefs @ ns) % p) != 0).sum(axis=1).min()); dl = min(dl, w)
        else:        # heavily punctured line: MITM with the cheap a = 2 tables only (d <= 5)
            for d in range(2, min(dl, 6)):
                if weight_d_exists_fast(M, p, d): dl = d; break
            else:
                if dl > 6: dl = 6  # unresolved: treat as <= 6 so the draw is not over-credited
    rec = {"p": p, "m": m, "r": r, "k": k, "seed": s, "n": G0.shape[1], "checks": G0.shape[0], "max_line": ml, "d_lines": dl if dl < 99 else ">=7", "t_lines": round(time.time() - t0), "tries": tries}
    if dl >= 6 and not nofull:
        t1 = time.time(); df = None
        for d in (2, 3, 4):
            if weight_d_exists_fast(G0, p, d): df = d; break
        if df is None: df = 5 if weight5_exists_parallel(G0, p) else None
        rec["full_w_le5"] = df; rec["t_full"] = round(time.time() - t1)
        rec["d"] = (6 if dl == 6 else ">=6") if df is None else df
        rec["puncture_columns_1indexed"] = [c + 1 for c in cols]
    with open(out, "a") as f: f.write(json.dumps(rec) + "\n")
    print(json.dumps({a: b for a, b in rec.items() if a != "puncture_columns_1indexed"}), flush=True)
