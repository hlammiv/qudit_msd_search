import json, glob, collections, sys
T = collections.defaultdict(lambda: [0.] * 4); rows = 0
for f in glob.glob(f"slots_d{sys.argv[1]}_w*.json"):
    for sup, sl, t in json.load(open(f))["rows"]:
        rows += 1
        for i in range(4): T[tuple(sorted(sl))][i] += t[i]
print("rows", rows)
for k, v in sorted(T.items()): print(k, 'u_T %.3f u_S %.3f acc_T %.2f acc_S %.2f' % tuple(v))
print('TOTAL u_T %.4f u_S %.4f acc_T %.2f acc_S %.2f' % tuple(sum(v[i] for v in T.values()) for i in range(4)))
