"""Pin the EXACT minimum distance of [[237,52,>=6]]_17 via the rebalanced MITM.

Run on lenore (125 GB / 32 cores):
    PYTHONPATH=$HOME/cap_validate_work venv/bin/python pin_d6.py
"""
import sys
import time

import numpy as np

from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.mindist_balanced import min_dependent_columns_balanced

PUNC = [6, 19, 32, 37, 40, 44, 47, 54, 59, 61, 69, 79, 84, 87, 95, 97, 103, 106, 109, 125,
        131, 132, 137, 138, 145, 151, 164, 177, 179, 185, 186, 187, 202, 204, 205, 206, 207,
        211, 212, 221, 246, 247, 248, 251, 254, 259, 266, 268, 269, 273, 283, 285]


def main():
    t0 = time.time()
    built = build_triorthogonal_code(17, 2, 10, PUNC)
    G0 = np.asarray(built["X_stab"]) % 17
    print(f"built [[{built['n']},{built['k']}]] full_rank={built['full_rank']} "
          f"dim(G0)={G0.shape[0]} shape={G0.shape}", flush=True)
    assert (built["n"], built["k"], built["full_rank"], G0.shape[0]) == (237, 52, True, 14)

    d = min_dependent_columns_balanced(G0, 17, d_max=6, n_jobs=32, ram_fraction=0.85)
    dt = time.time() - t0
    print(f"\nRESULT: min_distance([[237,52]]_17) = {d}", flush=True)
    print(f"VERDICT: d == 6 -> {'EXACTLY 6' if d == 6 else f'NOT 6 (={d})'}", flush=True)
    print(f"elapsed {dt/60:.1f} min", flush=True)
    sys.stdout.flush()


if __name__ == "__main__":
    main()
