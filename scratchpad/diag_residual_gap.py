"""Residual-frame failures of the closed-window escape vs their gaps (dump from escape_keep_rules.py with EKR_DUMP=1).
For kept (both bases) residual shots: failure value e = b - a mod d and its 2-adic valuation; the frame's X / Z parts
(valuations, weights); closed (gap2) and readout (gap) gaps.
    python -I scratchpad/diag_residual_gap.py dump.npz"""
import sys
from collections import Counter
import numpy as np
z = np.load(sys.argv[1]); d = int(z["d"])
def nu(v):
    v = np.asarray(v) % d; out = np.full(v.shape, np.log2(d).astype(int)); 
    for k in range(int(np.log2(d)) - 1, -1, -1):
        out[(v % (2 ** (k + 1)) != 0) & (v % (2 ** k) == 0)] = k
    return out
keep = z["keep_z"] & z["keep_x"]
e = (z["b"].astype(int) - z["a"].astype(int)) % d
g2 = np.minimum(z["gap2_z"], z["gap2_x"]); g = np.minimum(z["gap_z"], z["gap_x"])
F = z["frames"].astype(int) % d
print(f"d={d} residual shots {len(e)} kept {keep.sum()} fails {(keep & (e != 0)).sum()}")
for lo in (0, 4, 6, 7, 8, 9, 10):
    m = keep & (e != 0) & (g2 >= lo)
    print(f"  closed gap >= {lo:2d}: fails {m.sum():4d}  by e {dict(sorted(Counter(e[m]).items()))}  by nu(e) {dict(sorted(Counter(nu(e[m])).items()))}")
m = keep & (e != 0) & (g2 >= 6)
print("failing kept residual shots with closed gap >= 6:")
for i in np.nonzero(m)[0][:40]:
    xs, zs = F[i, :, 0], F[i, :, 1]
    print(f"  e={e[i]} gapZ/X closed {z['gap2_z'][i]:5.1f}/{z['gap2_x'][i]:5.1f} readout {z['gap_z'][i]:5.1f}/{z['gap_x'][i]:5.1f}"
          f"  X part {[int(v) for v in xs[xs != 0]]} @ {list(np.nonzero(xs)[0])}  Z part {[int(v) for v in zs[zs != 0]]} @ {list(np.nonzero(zs)[0])}")
# frame-type census among all kept residuals at closed gap >= 7: failing fraction by minimal valuation of the frame
mv = np.array([min([nu(v) for v in F[i].ravel() if v] or [np.log2(d)]) for i in range(len(F))])
for lo in (0, 7, 9):
    mk = keep & (g2 >= lo)
    print(f"  gap >= {lo}: kept by frame min-valuation {dict(sorted(Counter(mv[mk]).items()))}; fails {dict(sorted(Counter(mv[mk & (e != 0)]).items()))}")
