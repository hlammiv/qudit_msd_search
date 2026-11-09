"""Residual-frame priors for the escape decoder from an EKR_DUMP file: the K most frequent distinct nontrivial residual frames,
probability = count / accepted / p (unit location probability, as every DEM mechanism).
    python -I scratchpad/make_priors.py dump.npz accepted p K out.npz"""
import sys
import numpy as np
z = np.load(sys.argv[1]); acc, p, K = int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
F = z["frames"].astype(np.int64) % int(z["d"])
u, cnt = np.unique(F.reshape(len(F), -1), axis=0, return_counts=True)
o = np.argsort(-cnt)[:K]
np.savez(sys.argv[5], frames=u[o].reshape(-1, F.shape[1], 2).astype(np.int8), probs=cnt[o] / acc / p, counts=cnt[o])
print(len(u), "distinct frames; kept", len(o), "covering", cnt[o].sum(), "of", len(F), "; top counts", cnt[o][:8].tolist(), "min", cnt[o][-1])
