"""Exact distance of a level-set code (LEVEL_SETS.md) from coplanar parallel lines inside O.
The lower bound is d >= (p - b')p^(m-a'-1) from RM(m(p-1)-r+delta-1, m); for m = 4 it is (p - b') lines'
worth of points on a plane.  A minimum-weight word of the dual RM(m(p-1)-r-1, m) is (p - b) parallel lines
in a plane; if O contains delta = (b' - b) of those lines entirely, the word restricted off O has weight
(p - b - delta) p = the lower bound, so d is exact and the word is an explicit witness.
Search: for every direction v, count O-points on each line (fibres of the quotient by v); full lines are
points of the quotient F^4/<v> ~ F^3; delta full lines are coplanar (with direction v) iff their quotient
points are collinear.  usage: levelset_lines.py CODE.json OUT.json"""
import sys, json, itertools, time, numpy as np
from qmsd.reedmuller import points
src = json.load(open(sys.argv[1])); b = src.get("best", src); p, m, r = b["p"], b.get("m", 4), b["r"]; delta = b.get("delta", 3)
t0 = time.time()
P = np.asarray(points(m, p), dtype=np.int64); O = np.array(sorted(c - 1 for c in b["puncture_columns_1indexed"]))
inO = np.zeros(len(P), dtype=bool); inO[O] = True
a1, b1 = divmod(m * (p - 1) - r - 1, p - 1); a2, b2 = divmod(m * (p - 1) - r + delta - 1, p - 1)
assert a1 == a2 == 2 and m == 4, "this search handles the planar (a = 2) case"
nl_word, nl_bound = p - b1, p - b2
print(f"p={p} r={r} delta={delta}: dual min word = {nl_word} parallel lines in a plane; lower bound = {nl_bound} lines = {nl_bound*p}; need {nl_word-nl_bound} full lines of O", flush=True)
need = nl_word - nl_bound
dirs = [np.array(v) for v in itertools.product(range(p), repeat=m) if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1]
found = None; nfull_total = 0
for di, v in enumerate(dirs):
    # quotient coordinates: a basis of the complement; project x -> x - x_j v where j = first nonzero index of v, drop coord j
    j = int(np.nonzero(v)[0][0]); Q = (P - np.outer(P[:, j], v)) % p; Q = np.delete(Q, j, axis=1)
    key = Q @ (p ** np.arange(m - 1))
    cnt = np.bincount(key[inO], minlength=p ** (m - 1))
    full = np.nonzero(cnt == p)[0]; nfull_total += len(full)
    if len(full) >= need:
        pts3 = np.stack([(full // p ** i) % p for i in range(m - 1)], axis=1)
        # find `need` collinear points among the full-line quotient points
        S = {tuple(x) for x in pts3}
        for i1 in range(len(pts3)):
            for i2 in range(i1 + 1, len(pts3)):
                d3 = (pts3[i2] - pts3[i1]) % p
                line = [tuple((pts3[i1] + t * d3) % p) for t in range(p)]
                on = [q for q in line if q in S]
                if len(on) >= need:
                    found = {"direction": v.tolist(), "quotient_points": [list(map(int, q)) for q in on[:need]], "quotient_line_dir": d3.tolist(), "quotient_base": pts3[i1].tolist()}
                    break
            if found: break
    if found: break
    if di % 2000 == 0: print(f"  {di}/{len(dirs)} directions, {nfull_total} full lines so far ({time.time()-t0:.0f}s)", flush=True)
res = {"p": p, "r": r, "delta": delta, "lower_bound": nl_bound * p, "found": found, "full_lines_seen": int(nfull_total), "directions_scanned": di + 1, "seconds": round(time.time() - t0)}
if found:
    res["d_exact"] = nl_bound * p; res["note"] = f"{need} coplanar parallel lines of O -> explicit kernel word of weight {nl_bound*p} = lower bound"
json.dump(res, open(sys.argv[2], "w"), indent=1); print(json.dumps({k: v for k, v in res.items() if k != "found"}), "found" if found else "NOT FOUND", flush=True)
