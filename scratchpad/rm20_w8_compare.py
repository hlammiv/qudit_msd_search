"""Compare the span-3 dependent 8-sets of rm20_w8_count.py (residual method, recomputed here) with the union
of the brute-force shards (rm20_w8_brute.py).  usage: PYTHONPATH=. python scratchpad/rm20_w8_compare.py DIR [DIR ...]"""
import sys, glob, json
import numpy as np
src = open("scratchpad/rm20_w8_count.py").read()
head = src.split("for j in (2, 3, 4):")[0].replace('state = json.load(open(OUT)) if os.path.exists(OUT) else ', 'state = ')
g = {}; exec(compile(head, "rm20_w8_count.py", "exec"), g)
np_, p, m, W = g["np"], g["p"], g["m"], g["W"]
frame, flat = g["frame_and_flat"](m, 3, p)
others = np.array([x for x in flat if x not in set(frame)], dtype=np.int64)
M = g["frame_reduce"](frame); resid = M[4:, :][:, others].T.copy() % p
V, I, J = g["pair_dirs"](resid); kp = g["direction_keys"](V, p); ks = g["direction_keys"](resid, p)
triples = set()
for a, c in g["collisions"](kp, ks):
    x, y = int(I[a]), int(J[a])
    if c not in (x, y): triples.add(tuple(sorted((x, y, c))))
dep = set()
for tr in triples:
    for d in range(len(others)):
        if d not in tr: dep.add(tuple(sorted(tr + (d,))))
for a, b in g["collisions"](kp, kp):
    if a >= b: continue
    s = {int(I[a]), int(J[a]), int(I[b]), int(J[b])}
    if len(s) == 4: dep.add(tuple(sorted(s)))
brute = set(); shards = set(); sets = 0
for d in sys.argv[1:]:
    for f in glob.glob(f"{d}/span3_w*.json"):
        j = json.load(open(f)); shards.add(j["wid"]); sets += j["sets"]
        brute |= {tuple(x) for x in j["dependent"]}
nw = 20
mine = {s for s in dep if s[0] % nw in shards}
print(f"shards {sorted(shards)} ({len(shards)}/{nw}), brute sets {sets}, residual-method dependent {len(dep)} (in these shards {len(mine)}), brute dependent {len(brute)}")
print("agree:", mine == brute, "| only residual:", len(mine - brute), "| only brute:", len(brute - mine))
