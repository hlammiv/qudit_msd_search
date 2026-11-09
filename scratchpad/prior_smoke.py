"""Smoke test of Escape(priors=...): the {5,8,9,12} frame's detector pair gets an edge; decode that frame alone (p small)."""
import sys
import numpy as np
from primepower_msd.graft_run import VARIANTS
from primepower_msd.graft_zd import Escape
d = int(sys.argv[1]); z = np.load(sys.argv[2])
for pri in (None, (z["frames"].astype(np.int64), z["probs"])):
    ex = Escape(d, "X", **VARIANTS["one+gap-tick1-closed"], priors=pri)
    f = np.zeros((d - 1, 19, 2), dtype=np.int64)
    for v in range(1, d):
        f[v - 1, [5, 8, 9, 12], 1] = [v, -v, -v, v]
    f %= d
    V, O, _ = ex.prior_mechanisms(f, np.ones(d - 1))
    j = [tuple(int(x) for x in np.nonzero(V[0, ex.sel])[0])]
    print("priors" if pri else "none", "edges", len(ex.dem.edges), len(ex.dem_syn.edges), "pair edge", j, ex.dem_syn.edges.get(j[0]) and ex.dem_syn.edges[j[0]]["classes"][:2])
    r = ex.run(1e-3, d - 1, seed=1, decoder="peel-gap", inject=f, per_shot=True, chunk_seeded=True, flag_cert=True)
    print("   true obs", O.tolist(), "decoded flip", r["shot_flip"].tolist(), "closed gap", np.round(r["shot_gap_syn"], 1).tolist())
