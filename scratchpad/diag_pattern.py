"""Census of the Z-only residual frames supported on colour positions {5, 8, 9, 12} (escape dump): count, kept, failed, gaps.
    python -I scratchpad/diag_pattern.py dump.npz"""
import sys
import numpy as np
z = np.load(sys.argv[1]); d = int(z["d"]); F = z["frames"].astype(int) % d
S = [5, 8, 9, 12]
xs, zs = F[:, :, 0], F[:, :, 1]
supp = (zs != 0)
pat = ~np.any(xs != 0, axis=1) & np.all(supp[:, S], axis=1) & (supp.sum(axis=1) == 4)
keep = z["keep_z"] & z["keep_x"]; e = (z["b"].astype(int) - z["a"].astype(int)) % d
g2x = z["gap2_x"]; g2z = z["gap2_z"]
print(f"d={d} residual frames {len(F)}; pattern frames {pat.sum()}; kept {(pat & keep).sum()}; failed kept {(pat & keep & (e != 0)).sum()}")
v = np.unique(np.round(g2x[pat & keep], 1), return_counts=True)
print("  closed X gap of kept pattern frames:", dict(zip(v[0].tolist(), v[1].tolist())))
print("  failed kept with min gap >= 7:", int((pat & keep & (e != 0) & (np.minimum(g2x, g2z) >= 7)).sum()),
      "of all failed kept >= 7:", int((keep & (e != 0) & (np.minimum(g2x, g2z) >= 7)).sum()))
print("  pattern values (first 10):", [tuple(zs[i, S]) for i in np.nonzero(pat)[0][:10]])
