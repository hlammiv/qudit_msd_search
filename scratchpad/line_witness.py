"""Find an explicit weight-6 line-supported word of an m = 2 line-spread code from a certified JSONL
record (puncture_columns_1indexed), searching the most-punctured lines first (exact MITM on the line).
Prints the witness or reports that no line carries a weight-6 word (then d >= 6 only).
usage: line_witness.py RECORD.jsonl"""
import sys, json, itertools, numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.structured_ad import _rank, mitm_streamed
from qmsd.triorthogonal import build_triorthogonal_code
r = [json.loads(l) for l in open(sys.argv[1]) if l.strip()][0]
p, rr, cols = r["p"], r["r"], r["puncture_columns_1indexed"]
P = np.asarray(points(2, p), dtype=np.int64); G = np.asarray(rm_generator(rr, 2, p), dtype=np.int64) % p
G0 = np.asarray(build_triorthogonal_code(p, 2, rr, cols, G=G)["G0"], dtype=np.int64) % p
cs = set(c - 1 for c in cols); keep = [l for l in range(len(P)) if l not in cs]; pos = {l: j for j, l in enumerate(keep)}
lines = []
for a in range(p):
    key = (P[:, 0] + a * P[:, 1]) % p; lines += [np.flatnonzero(key == c) for c in range(p)]
lines += [np.flatnonzero(P[:, 1] == c) for c in range(p)]
jc = [sum(1 for i in ln if i in cs) for ln in lines]
print(f"p={p} k={len(cols)} n={G0.shape[1]} checks={G0.shape[0]} max outputs on a line {max(jc)}", flush=True)
for ln, j in sorted(zip(lines, jc), key=lambda x: -x[1]):
    cc = np.array([pos[i] for i in ln if i in pos], dtype=np.int64)
    if len(cc) < 6: continue
    M = G0[:, cc]; M = M[[t for t in range(M.shape[0]) if M[t].any()]]
    if _rank(M.T.copy(), p) == len(cc): continue
    for d in range(2, 7):
        ws = list(itertools.islice(mitm_streamed(M, d, p), 1))
        if ws:
            supp, vals = ws[0]; g0 = [int(cc[s]) for s in supp]; coef = [int(v) for v in vals]
            vec = np.zeros(G0.shape[1], dtype=np.int64); vec[g0] = coef; assert ((G0 @ vec) % p == 0).all()
            print(json.dumps({"witness_weight": d, "witness_G0_cols": g0, "witness_coeffs": coef, "outputs_on_line": j}), flush=True)
            sys.exit(0)
print("NO line word of weight <= 6 found on any line: d >= 6 only from this check", flush=True)
