import sys, json, itertools as it
d = int(sys.argv[1]); K = int(sys.argv[2]); NW = int(sys.argv[3]); sys.argv = ['x', str(d)]
exec(open('proxy_validate.py').read().split('out = {')[0])
sups = set()
for sup in it.combinations(range(n), 3):
    for Hc in (MZ, H):
        for v in it.product(range(1, d), repeat=3):
            e = np.zeros(n, dtype=np.int64); e[list(sup)] = v
            if not ((Hc @ e) % d).any(): sups.add(sup)
sups = sorted(sups)
def task(sup, sl):
    t = [0., 0., 0., 0.]
    for pv in it.product(paulis, repeat=3):
        fs = [(s, i, a, b) for s, i, (a, b) in zip(sl, sup, pv)]
        aT, uT = run("T", fault_set(fs)); aS, uS = run("S", fault_set(fs))
        t[0] += uT; t[1] += uS; t[2] += aT; t[3] += aS
    return [list(sup), list(sl), t]
import time
t0 = time.time()
allt = [(sup, sl) for sup in sups for sl in it.product(range(4), repeat=3)]
res = []
for j, a in enumerate(allt):
    if j % NW == K:
        res.append(task(*a))
        json.dump({"d": d, "n_supports": len(sups), "rows": res}, open(f"slots_d{d}_w{K}.json", "w"))
print("worker", K, len(res), round(time.time() - t0), flush=True)
