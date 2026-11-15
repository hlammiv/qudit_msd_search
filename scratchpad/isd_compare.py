"""Compare the union of ISD d''-words (isd_paired_w8.py) with the structural line-pair family of census_paired_w8.py.
usage: PYTHONPATH=. python scratchpad/isd_compare.py SWEEP.json T SAMPLE ISD1.json [ISD2.json ...]"""
import sys, json, glob, itertools, numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.triorthogonal import shorten_matrix
from qmsd.structured_ad import _rref, _null_space
src, t, s = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); p = 7
Y = [d for d in json.load(open(src))["draws"] if d["t"] == t and d["sample"] == s][0]["Y"]; Yset = set(Y)
P3 = np.asarray(points(3, p), dtype=np.int64); idx3 = {tuple(int(v) for v in x): i for i, x in enumerate(P3)}
offY = [x for x in range(p ** 3) if x not in Yset]; posY = {x: j for j, x in enumerate(offY)}
R6 = np.asarray(rm_generator(6, 3, p), dtype=np.int64) % p; R, pv = _rref(np.asarray(shorten_matrix(R6, Y, p)) % p, p); H6 = R[:len(pv)]
_dirs = sorted({min(tuple(int(v) for v in (c * np.array(d)) % p) for c in range(1, p)) for d in itertools.product(range(p), repeat=3) if any(d)})
LINES = sorted({tuple(sorted(idx3[tuple(int(v) for v in (P3[a] + c * np.array(d)) % p)] for c in range(p))) for a in range(p ** 3) for d in _dirs})
cand = [L for L in LINES if sum(x in Yset for x in L) >= 2]
struct = set()
for A, B in itertools.combinations(cand, 2):
    off = sorted(x for x in set(A) | set(B) if x not in Yset)
    if len(off) != 8: continue
    for v in np.asarray(_null_space(H6[:, [posY[x] for x in off]], p)) % p:
        if np.count_nonzero(v) == 8:
            inv = pow(int(v[0]), p - 2, p); struct.add(tuple(zip(off, (int(x) * inv % p for x in v))))
isd = {}
for f in sys.argv[4:]:
    for w in json.load(open(f))["words"]:
        isd.setdefault(tuple(tuple(c) for c in w), []).append(f)
U = set(isd)
print(f"t={t}: structural {len(struct)}, ISD union {len(U)} over {len(sys.argv) - 4} seeds, ISD-only {len(U - struct)}, structural-only {len(struct - U)}")
per = [sum(1 for w in struct if f in isd.get(w, [])) for f in sys.argv[4:]]
print("structural words found per seed:", per)
miss = struct - U
for w in list(U - struct)[:3]: print("ISD-only example:", w)
json.dump({"t": t, "structural": len(struct), "isd_union": len(U), "isd_only": [list(map(list, w)) for w in U - struct],
           "structural_only": len(miss)}, open(f"runs/weighted/isd/compare_t{t}.json", "w"), indent=1)
