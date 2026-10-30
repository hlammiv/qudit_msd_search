"""Quick diagnostic: 2+3 probe cost vs left-table size (cache-resident -> DRAM)."""
import math
import time

import numpy as np

from qmsd import mindist_nb as m

p, r = 7, 31
for n, cols in ((60, (0, 60)), (100, (0, 60)), (400, (0, 12))):
    rng = np.random.default_rng(1)
    H = rng.integers(0, p, size=(r, n))
    pw = m._powers_u64(p, r)
    lc, ls = m._build_left(H, p, 2, pw)
    bs, sh = m._left_index(lc)
    m._scan_right_b3_idx(H, p, pw, lc, ls, bs, sh, 0, 1)
    t = time.perf_counter()
    m._scan_right_b3_idx(H, p, pw, lc, ls, bs, sh, *cols)
    dt = time.perf_counter() - t
    npb = sum(math.comb(n - 1 - i, 2) for i in range(*cols)) * 36
    print(f"n={n} N={lc.shape[0]} table={(lc.nbytes + ls.nbytes) >> 20} MB probes={npb:.3g} "
          f"{npb / dt:.3g}/s  {1e9 * dt / npb:.0f} ns/probe")
