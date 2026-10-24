import sys, json, time, itertools as it, numpy as np
from primepower_msd.meas_distill import signed_rm
from primepower_msd.meas_sim import Protocol, count_input_failures_fast, count_inner_failures, min_weight_logicals
d = int(sys.argv[1]); out = {}
for m in range(3, 8):
    S, tau = signed_rm(1, m); n = S.shape[1]; k = n - 2 * S.shape[0]
    bc = np.array(list(it.combinations(range(k + 1), 2)), dtype=np.int64)
    P = Protocol(S, tau, d, bc); t = time.time()
    for w in (1, 2):
        assert not min_weight_logicals(P, w)[1], (m, w)
        assert not count_input_failures_fast(P, w, jobs=4), (m, w)
    inner = count_inner_failures(P, 3); inp = count_input_failures_fast(P, 3, jobs=4)
    C3 = inner + sum(inp.values())
    rec = dict(n=n, k=k, inputs=P.n_in, checks=P.M, L=P.locations(), T_per_out=1 + 2 * P.s * n / k, C3=C3,
               input={str(a): b for a, b in inp.items()}, inner=inner, coef=C3 / P.n_in / (d - 1) ** 3)
    out[str(m)] = rec; print(d, m, json.dumps(rec), f"{time.time()-t:.0f}s", flush=True)
json.dump(out, open(f"runs/meas_distill/order3_d{d}.json", "w"), indent=1)
