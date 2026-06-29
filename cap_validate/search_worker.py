import sys, json
from qmsd.search import random_search
# args: points("m:k,m:k,..."), outfile, trials, seed
points = [tuple(int(x) for x in p.split(":")) for p in sys.argv[1].split(",")]
outfile, trials, seed = sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
out = {}
for m, k in points:
    N = 3**m
    codes = random_search(3, m, trials=trials, seed=seed, target_k=k, sampler="uniform", n_jobs=-1, max_distance=5)
    codes += random_search(3, m, trials=trials//3, seed=seed+1, target_k=k, sampler="capset_climb", n_jobs=-1, max_distance=5)
    pool = [x for x in codes if x.d and x.d >= 2]
    if not pool:
        print(f"m={m} k={k}: none", flush=True); continue
    dmax = max(x.d for x in pool)
    punc = []
    seen = set()
    for x in pool:
        if x.d == dmax and x.puncture_columns:
            t = tuple(x.puncture_columns)
            if t not in seen:
                seen.add(t); punc.append(list(t))
    out[f"{m}:{k}"] = {"m": m, "k": k, "n": N - k, "d": dmax, "n_maxd": len(punc), "maxd_punctures": punc[:25]}
    print(f"m={m} k={k}: d_max={dmax}  ({len(punc)} distinct max-d codes)", flush=True)
json.dump(out, open(outfile, "w"))
print("WROTE", outfile, flush=True)
