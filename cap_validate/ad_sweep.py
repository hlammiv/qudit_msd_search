import sys, json
from joblib import Parallel, delayed
from qmsd.search import random_search
from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.weightdist import exact_distance_and_Ad
from qmsd.reedmuller import rm_generator, r_max

# args: points("m:k:d,..."), trials, max_ad, ad_jobs, outfile, seed
points = [tuple(int(x) for x in p.split(":")) for p in sys.argv[1].split(",")]
trials, max_ad, ad_jobs, outfile, seed = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], int(sys.argv[6])
PAPER = {(4,9,3):648,(5,43,3):1700,(5,37,4):880,(6,62,4):3972,(5,28,5):1104}

def ad_of(punc, m, k, d):
    r = r_max(m,3); G = rm_generator(r,m,3)
    b = build_triorthogonal_code(3,m,r,punc,G=G)
    res = exact_distance_and_Ad(b["X_stab"], 3, max_words=30_000_000)
    return res.get("B_d") if res.get("distance")==d else None

results = {}
for m,k,d in points:
    N = 3**m
    cap_dist = min(d+1, 6)
    codes  = random_search(3,m,trials=trials,   seed=seed,   target_k=k, sampler="capset_climb", n_jobs=-1, max_distance=cap_dist)
    codes += random_search(3,m,trials=trials//2, seed=seed+1, target_k=k, sampler="uniform",      n_jobs=-1, max_distance=cap_dist)
    seen, pool = set(), []
    for x in codes:
        if x.d==d and x.puncture_columns:
            t = tuple(x.puncture_columns)
            if t not in seen: seen.add(t); pool.append(list(t))
    ads = Parallel(n_jobs=ad_jobs)(delayed(ad_of)(p,m,k,d) for p in pool[:max_ad])
    valid = [(a,p) for a,p in zip(ads, pool[:max_ad]) if a]
    min_ad, best = (min(valid, key=lambda x:x[0]) if valid else (None,None))
    paper = PAPER.get((m,k,d))
    imp = round(paper/min_ad,2) if (paper and min_ad) else None
    results[f"{m}:{k}:{d}"] = {"m":m,"k":k,"d":d,"n":N-k,"n_codes":len(pool),"checked":len(valid),
                               "min_A_d":min_ad,"paper_A_d":paper,"improvement":imp,"best_punc":best}
    json.dump(results, open(outfile,"w"))
    print(f"[[{N-k},{k},{d}]]_3: {len(pool)} d={d} codes, min A_d={min_ad} vs paper {paper}  ({imp}x better)", flush=True)
print("DONE", flush=True)
