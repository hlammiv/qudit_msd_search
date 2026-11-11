"""Brute-force cross-check of rm20_w8_count.py: batched F_7 rank of the 35 x 8 cubic evaluation
matrices of frame sets.  span 3: every set {frame} + 4 of the 339 others, split by the first
extra index (x % NW == WID); span 4: SAMPLE random 3-subsets of the 2396 others.
Writes the dependent sets (indices into `others`, sorted) as JSON.

usage: PYTHONPATH=. python scratchpad/rm20_w8_brute.py SPAN WID NW OUT [SAMPLES]
"""
import sys, os, json, time, itertools
import numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.pointsets import batched_rank_modp, frame_and_flat

span, wid, nw, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
p, m = 7, 4
R3 = np.asarray(rm_generator(3, m, p), dtype=np.int64) % p
frame, flat = frame_and_flat(m, span, p)
others = np.array([x for x in flat if x not in set(frame)], dtype=np.int64)
extra = 8 - (span + 1)
E = R3.T                                                    # (2401, 35)
BATCH = int(os.environ.get('BATCH', '50000'))
T0 = time.time(); dep = []; nsets = 0
def test(ch):
    global nsets
    T = np.concatenate([np.broadcast_to(np.array(frame), (len(ch), span + 1)), others[ch]], axis=1)
    rk = batched_rank_modp(E[T], p)
    nsets += len(ch)
    for row in ch[rk < 8]: dep.append(sorted(int(v) for v in row))
if span == 3:
    n = len(others)
    for x in range(wid, n, nw):
        it = itertools.combinations(range(x + 1, n), extra - 1)
        while True:
            r = np.fromiter(itertools.chain.from_iterable(itertools.islice(it, BATCH)), dtype=np.int64)
            if r.size == 0: break
            r = r.reshape(-1, extra - 1)
            test(np.concatenate([np.full((len(r), 1), x), r], axis=1))
        print(f"[{time.time()-T0:6.0f}s] x={x} sets {nsets} dependent {len(dep)}", flush=True)
else:
    samples = int(sys.argv[5]); rng = np.random.default_rng(1000 + wid); n = len(others)
    while nsets < samples:
        ch = np.sort(np.array([rng.choice(n, extra, replace=False) for _ in range(BATCH)]), axis=1)
        test(ch)
        if nsets % 1_000_000 == 0: print(f"[{time.time()-T0:6.0f}s] sets {nsets} dependent {len(dep)}", flush=True)
json.dump({"span": span, "wid": wid, "nw": nw, "sets": nsets, "dependent": dep, "seconds": round(time.time() - T0)}, open(out, "w"))
print("done", nsets, len(dep), flush=True)
