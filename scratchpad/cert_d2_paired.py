"""d''(Y) >= 8 test for a paired-cap draw (probe_paired_cap.py): H = RM(6,3) shortened at Y (343 - t columns);
weights <= 5 (serial + parallel 2+3), 6 (3+3 witness) and 7 (3+4 witness).  With d'(Y) = 4 this makes the
[[686 - 2t, 2t]] single-species code d = min(2 d', d'') = 8.
Checkpointed: OUT.json.ckpt records the finished stages and the weight-7 blocks done (keyed by a fingerprint of H);
rerunning the same command resumes from it.
usage: PYTHONPATH=. python scratchpad/cert_d2_paired.py SWEEP.json T SAMPLE OUT.json"""
import sys, os, json, time, hashlib, numpy as np
from qmsd.reedmuller import rm_generator
from qmsd.triorthogonal import shorten_matrix
from qmsd.mindist_nb import weight_d_exists_fast, weight5_exists_parallel, weight6_witness_3x3, weight7_witness_3x4
src, t, s, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
p = 7; BLOCK = 8; ck_path = out + ".ckpt"
dr = [d for d in json.load(open(src))["draws"] if d["t"] == t and d["sample"] == s][0]
R6 = np.asarray(rm_generator(6, 3, p), dtype=np.int64) % p
H = np.asarray(shorten_matrix(R6, dr["Y"], p), dtype=np.int64) % p
fp = hashlib.sha256(np.ascontiguousarray(H).tobytes()).hexdigest()
ck = {"fingerprint": fp, "low_done": False, "d2_low": None, "w6_done": False, "w6_witness": None, "w7_next": 0, "w7_seconds": 0.0}
if os.path.exists(ck_path):
    old = json.load(open(ck_path))
    if old.get("fingerprint") == fp:
        ck = old; print(f"resuming from {ck_path}: {ck}", flush=True)
    else:
        print(f"ignoring {ck_path}: different matrix", flush=True)

def save():
    tmp = ck_path + ".tmp"; json.dump(ck, open(tmp, "w")); os.replace(tmp, ck_path)

res = {"t": t, "sample": s, "Y": dr["Y"], "d1": dr["d1"], "cols": int(H.shape[1])}
T0 = time.time()
if not ck["low_done"]:
    d2 = None
    for w in (2, 3, 4):
        if weight_d_exists_fast(H, p, w): d2 = w; break
    if d2 is None and weight5_exists_parallel(H, p): d2 = 5
    ck.update(low_done=True, d2_low=d2); save()
d2 = ck["d2_low"]
if d2 is None and not ck["w6_done"]:
    w6 = weight6_witness_3x3(H, p)
    ck.update(w6_done=True, w6_witness=w6); save()
res["w6_witness"] = ck["w6_witness"]
if d2 is None and ck["w6_witness"] is not None: d2 = 6
print(f"weights <= 6 done: d''={d2} [{time.time()-T0:.0f}s]", flush=True)
if d2 is None:
    def log(a, b, sec):
        ck["w7_next"] = b; ck["w7_seconds"] += sec; save()
        print(f"  w7 lead cols {a}-{b} {sec:.0f}s", flush=True)
    w7 = weight7_witness_3x4(H, p, block=BLOCK, log=log, i_min=ck["w7_next"])
    res["w7_witness"] = w7
    if w7 is not None: d2 = 7
res["d2"] = d2 if d2 is not None else ">=8"; res["seconds"] = time.time() - T0
res["d_code"] = min(2 * dr["d1"], d2) if d2 is not None else min(2 * dr["d1"], 8)
print({k: v for k, v in res.items() if k != "Y"}, flush=True)
json.dump(res, open(out, "w"), indent=1)
