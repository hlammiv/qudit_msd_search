"""Validate the COMPILED dwd_core (C+OpenMP) against the trusted qmsd engine.

Run from /home/hlamm/Desktop/QC/prime_msd:
    python -m dist_weightdist.validate_core

PROTOCOL (design sec.9, extending validate.py to the real compiled core path):
  For each oracle code [[206,37,4]]_3 (K=14) and [[667,62,4]]_3 (K=16):
    1. build G0 = build_triorthogonal_code(...)["X_stab"]; write it as .g0.
    2. run dwd_core to histogram all 3^K codewords (exercising the DISTRIBUTED path:
       partition into 3^t blocks, run as MULTIPLE chunks/partials, merge them).
    3. assert merged A == qmsd.weightdist.weight_enumerator EXACTLY (vector equality).
    4. run merge.certify -> assert d and A_d match qmsd.exact_distance_and_Ad AND the
       published oracle values (A_d = 880, 3972); assert all dual invariants.
    5. assert dwd_core's own packed-vs-scalar selfcheck passes.

An exact match on merged A guarantees an exact match on (d, A_d) since the engine
feeds the SAME extract_d_and_Ad path validated in validate.py.
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np

from qmsd import weightdist as qwd
from qmsd.oracle import load_oracle
from qmsd.triorthogonal import build_triorthogonal_code

from . import correctness as cx
from . import merge as mg
from . import partition as pt
from .harness import core_binary, run_chunk
from .io_g0 import write_g0

ORACLES = [("[[206,37,4]]_3", 4, 880), ("[[667,62,4]]_3", 4, 3972)]


def build_g0(label: str):
    oc = next(o for o in load_oracle() if o.label == label)
    built = build_triorthogonal_code(oc.p, oc.m, oc.r_max, oc.puncture_columns_1indexed)
    G0 = np.asarray(built["X_stab"])
    return G0, oc.p


def selfcheck(core: Path) -> None:
    res = subprocess.run([str(core), "selfcheck"], capture_output=True, text=True)
    assert res.returncode == 0, f"dwd_core selfcheck failed:\n{res.stdout}\n{res.stderr}"
    assert "SELFCHECK PASSED" in res.stdout, res.stdout
    print("[ok] dwd_core packed-vs-scalar selfcheck PASSED")


def check_oracle(core: Path, label: str, exp_d: int, exp_Ad: int) -> None:
    G0, q = build_g0(label)
    assert q == 3
    K, n = G0.shape

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        g0path = td / "G0.g0"
        write_g0(g0path, G0, q)

        # DISTRIBUTED path: split the 3^K space into 3^t blocks and run as several
        # chunks -> several .partial files (exactly what the 2-node job produces).
        t = min(2, K)                       # 9 blocks -> 3 chunks of 3 blocks
        nblk = pt.pow3(t)
        chunks = pt.chunkify(0, nblk, max(1, nblk // 3))
        partials = []
        for i, (s, c) in enumerate(chunks):
            pp = td / f"chunk{i}.partial"
            run_chunk(core, str(g0path), t, s, c, str(pp))
            partials.append(str(pp))

        # merged A from partials (the SAME merge the distributed engine uses)
        A, nn, KK = mg.merge_to_A(partials)
        assert nn == n and KK == K

    # 3. exact vector match against qmsd's independent enumeration (the one slow ref enum)
    A_ref = qwd.weight_enumerator(G0, q, max_words=q ** K + 1)
    assert A == list(A_ref), f"{label}: dwd_core A != qmsd.weight_enumerator"

    # 4. certify (d, A_d, B) and compare to qmsd's MacWilliams on the SAME (validated) A_ref
    #    -- this reproduces qmsd.exact_distance_and_Ad exactly without a 2nd 3^K enumeration.
    res = cx.extract_d_and_Ad(A, q)
    cx.assert_dual_invariants(res["B"], A, q)
    ref_B = qwd.macwilliams(list(A_ref), q)
    ref_d = next(w for w in range(1, len(ref_B)) if ref_B[w] > 0)
    assert res["d"] == ref_d == exp_d, (label, res["d"], ref_d, exp_d)
    assert res["A_d"] == ref_B[ref_d] == exp_Ad, (label, res["A_d"], ref_B[ref_d], exp_Ad)
    assert res["B"] == ref_B, f"{label}: dual B disagrees with qmsd.macwilliams"
    assert res["K"] == K
    print(f"[ok] {label}: dwd_core A == qmsd EXACTLY; K={K} n={n} d={res['d']} "
          f"A_d={res['A_d']} (== published {exp_Ad}); invariants pass "
          f"(B_0=1, sum_w B_w = 3^{n - K}, all B_w>=0)")


def main():
    core = core_binary()
    print(f"[core] {core}")
    selfcheck(core)
    for label, d, Ad in ORACLES:
        check_oracle(core, label, d, Ad)
    print("ALL dwd_core VALIDATION CHECKS PASSED (compiled core == qmsd EXACTLY)")


if __name__ == "__main__":
    main()
