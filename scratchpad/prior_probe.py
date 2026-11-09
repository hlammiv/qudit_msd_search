"""Detector signature of the {5,8,9,12} Z residual frame in the X-basis escape (closed window) and how the DEM sees it."""
import sys
import numpy as np
from primepower_msd.graft_run import VARIANTS
from primepower_msd.graft_zd import Escape
d = int(sys.argv[1])
ex = Escape(d, "X", **VARIANTS["one+gap-tick1-closed"])
f = np.zeros((1, 19, 2), dtype=np.int64); f[0, [5, 8, 9, 12], 1] = [1, -1, -1, 1]
V, O, _ = ex.prior_mechanisms(f % d, [1.0])
for name, sel in (("readout dem", ex.sel), ("closed dem", ex.sel_syn)):
    nz = [(int(j), int(V[0, sel][j]), ex.mem.det_tags[sel[j]]) for j in np.nonzero(V[0, sel])[0]]
    print(name, "fires", nz, "obs", int(O[0]))
    if len(nz) == 2:
        dem = ex.dem if name.startswith("readout") else ex.dem_syn
        e = dem.edges.get((nz[0][0], nz[1][0]))
        print("   existing edge:", e and {"pat": e["pat"], "classes": e["classes"][:3], "P": e["P"]})
print("masked (postselected) dets fired:", int((V[0, ex.mask] != 0).sum()))
