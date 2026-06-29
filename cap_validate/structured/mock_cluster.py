"""LOCAL two-mock-node test of the Method-2 distribution logic.

NEVER contacts a remote host.  Simulates the local(20-core)+lenore(32-core) split as two
LOCAL worker subprocesses over disjoint direction slices through a shared temp directory,
then merges and asserts the merged result equals the single-process whole-family result.
Also exercises checkpoint/resume by deleting one node's partials and re-running.

Run:  python -m cap_validate.structured.mock_cluster
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

from . import capcode as cc
from . import distribute as dist
from . import enum_families as ef


def run(family="w18"):
    fam = ef.FAMILIES[family]
    Spts = cc.load_cap()["cap_pts"]
    Hs = ef.directions(fam)
    n_dirs = Hs.shape[0]

    # ground truth: single process over the whole family
    truth = ef.family_min_punct(fam, Spts, Hs)
    print(f"[mock] family={family} n_dirs={n_dirs}  single-process "
          f"min_punct={truth['min_punct']} (max_inter={truth['max_inter']})")

    # split across two mock nodes weighted 20:32 (local cores : lenore cores), chunked
    plan = dist.split_blocks(family, n_dirs, node_weights=[20.0, 32.0], chunk=12000)
    dist.assert_cover(plan, n_dirs)
    node_names = ["local", "lenore_mock"]
    print(f"[mock] plan: " +
          ", ".join(f"{node_names[ni]}={sum(b.count for b in plan[ni])}dirs"
                    f"/{len(plan[ni])}blocks" for ni in plan))

    with tempfile.TemporaryDirectory() as td:
        # 1. each node runs its blocks in a SEPARATE subprocess (true process split)
        for ni, blocks in plan.items():
            argv = [sys.executable, "-m", "cap_validate.structured.distribute",
                    node_names[ni], family]
            for b in blocks:
                argv += [str(b.start), str(b.count)]
            argv += [td]
            subprocess.run(argv, check=True, cwd=str(cc.REPO))

        merged = dist.merge_partials(td, family, n_dirs)
        print(f"[mock] merged min_punct={merged['min_punct']} "
              f"max_inter={merged['max_inter']} witness_dir={merged['witness_dir']} "
              f"blocks={merged['n_blocks']}")
        assert merged["min_punct"] == truth["min_punct"]
        assert merged["max_inter"] == truth["max_inter"]

        # 2. checkpoint/resume: wipe lenore_mock's partials, re-run, re-merge -> identical
        for p in Path(td).glob("lenore_mock_*.partial"):
            p.unlink()
        argv = [sys.executable, "-m", "cap_validate.structured.distribute",
                "lenore_mock", family]
        for b in plan[1]:
            argv += [str(b.start), str(b.count)]
        argv += [td]
        subprocess.run(argv, check=True, cwd=str(cc.REPO))
        merged2 = dist.merge_partials(td, family, n_dirs)
        assert merged2["min_punct"] == truth["min_punct"]
        assert merged2["max_inter"] == truth["max_inter"]
        print("[mock] resume OK: merged identical after wiping+rerunning one node")

    print("[mock] PASS: distributed == single-process; coverage exact; resume idempotent")
    return {"family": family, "min_punct": merged["min_punct"], "truth": truth}


if __name__ == "__main__":
    run("w18")
    run("w27")
