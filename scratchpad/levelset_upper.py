"""Best planar upper bound for a level-set code when no exact coplanar-lines witness exists: over planes
through any full line of O (direction v, quotient point q, every quotient line through q), place the
(p - b) parallel lines of a minimum-weight dual word on the lines with the most O-points.  Kernel weight =
(p - b) p - (O-points covered).  Writes the best plane in levelset_lines.py's format, so levelset_witness.py
can build and verify the word.  usage: levelset_upper.py CODE.json OUT.json"""
import sys, json, itertools, time, numpy as np
from qmsd.reedmuller import points
src = json.load(open(sys.argv[1])); b = src.get("best", src); p, m, r = b["p"], 4, b["r"]; t0 = time.time()
P = np.asarray(points(m, p), dtype=np.int64); O = np.array(sorted(c - 1 for c in b["puncture_columns_1indexed"]))
inO = np.zeros(len(P), dtype=bool); inO[O] = True
a1, b1 = divmod(m * (p - 1) - r - 1, p - 1); assert a1 == 2; nl = p - b1
dirs = [np.array(v) for v in itertools.product(range(p), repeat=m) if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1]
d3s = [np.array(v) for v in itertools.product(range(p), repeat=3) if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1]
best = (-1, None); nfull = 0
for v in dirs:
    j = int(np.nonzero(v)[0][0]); Q = (P - np.outer(P[:, j], v)) % p; Q = np.delete(Q, j, axis=1)
    key = Q @ (p ** np.arange(3)); cnt = np.bincount(key[inO], minlength=p ** 3)
    full = np.nonzero(cnt == p)[0]; nfull += len(full)
    for fk in full:
        q = np.array([(fk // p ** i) % p for i in range(3)])
        for d3 in d3s:
            keys = [int((((q + t * d3) % p) @ (p ** np.arange(3)))) for t in range(p)]
            c = np.sort(cnt[keys])[::-1][:nl]; s = int(c.sum())
            if s > best[0]:
                full_on = [list(map(int, (q + t * d3) % p)) for t in range(p) if cnt[keys[t]] == p]
                best = (s, {"direction": v.tolist(), "quotient_base": q.tolist(), "quotient_line_dir": d3.tolist(), "quotient_points": full_on})
res = {"p": p, "r": r, "word_lines": nl, "full_lines_total": int(nfull), "max_O_points_on_word": best[0], "upper_bound": nl * p - best[0], "found": best[1],
       "lower_bound": b.get("d_lower_bound"), "seconds": round(time.time() - t0)}
json.dump(res, open(sys.argv[2], "w"), indent=1); print({k: v for k, v in res.items() if k != "found"}, flush=True)
