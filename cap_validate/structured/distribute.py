"""Distribute the Method-2 enumeration across blocks / nodes (and a local mock cluster).

Unit of work: a contiguous slice [start, start+count) of the direction list of a family
(the directions are independent; max_topj over a slice is a pure reduction).  Merge rule
is MAX of per-block max|supp cap S|  ->  MIN of punctured weight.  This commutes with any
partition, so splitting across 2 machines and merging gives the identical global result.

Provenance: each partial records the family, slice, block max, and the witness direction,
plus a checksum (count of directions visited) so a merge can assert full coverage.

NO host is contacted here.  ``mock_cluster`` runs two LOCAL worker processes over disjoint
slices through a shared directory; the real 2-machine launch is the RUNBOOK the user fills.
"""
from __future__ import annotations

import json

import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import capcode as cc
from . import enum_families as ef
from . import subspaces as sub


# --------------------------------------------------------------------------- plan
@dataclass(frozen=True)
class Block:
    family: str
    start: int
    count: int


def split_blocks(family: str, n_dirs: int, node_weights: list[float],
                 chunk: int) -> dict[int, list[Block]]:
    """Partition [0, n_dirs) across nodes (proportional to node_weights), each node's
    span further cut into <= ``chunk`` sized blocks.  Returns {node_idx: [Block,...]}.
    Guarantees an exact, disjoint cover of [0, n_dirs)."""
    w = np.asarray(node_weights, dtype=float)
    bounds = np.floor(np.cumsum(w) / w.sum() * n_dirs).astype(int)
    bounds[-1] = n_dirs
    plan: dict[int, list[Block]] = {}
    lo = 0
    for ni, hi in enumerate(bounds):
        blocks = []
        s = lo
        while s < hi:
            c = min(chunk, hi - s)
            blocks.append(Block(family, s, c))
            s += c
        plan[ni] = blocks
        lo = hi
    return plan


def assert_cover(plan: dict[int, list[Block]], n_dirs: int):
    seen = np.zeros(n_dirs, dtype=np.int8)
    for blocks in plan.values():
        for b in blocks:
            seen[b.start:b.start + b.count] += 1
    assert seen.min() == 1 and seen.max() == 1, "plan is not a disjoint exact cover"


# --------------------------------------------------------------- worker (per block)
def run_block(family: str, start: int, count: int, Spts=None, Hs=None) -> dict:
    fam = ef.FAMILIES[family]
    if Spts is None:
        Spts = cc.load_cap()["cap_pts"]
    if Hs is None:
        Hs = ef.directions(fam)
    return ef.family_block_max(fam, Spts, Hs, start, count)


def run_node_jobs(blocks: list[Block], out_dir: str, node_name: str) -> list[str]:
    """Run all blocks for one node, writing one .partial JSON each.  Resumable: a block
    whose .partial already exists (and validates) is skipped."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cap = cc.load_cap()
    Spts = cap["cap_pts"]
    # cache direction matrices per family
    hs_cache: dict[str, np.ndarray] = {}
    written = []
    for b in blocks:
        pf = out / f"{node_name}_{b.family}_{b.start}_{b.count}.partial"
        if pf.exists():
            try:
                d = json.loads(pf.read_text())
                if d.get("count") == b.count and d.get("start") == b.start:
                    written.append(str(pf))
                    continue
            except Exception:
                pass
        if b.family not in hs_cache:
            hs_cache[b.family] = ef.directions(ef.FAMILIES[b.family])
        res = run_block(b.family, b.start, b.count, Spts, hs_cache[b.family])
        res["node"] = node_name
        pf.write_text(json.dumps(res))
        written.append(str(pf))
    return written


# --------------------------------------------------------------------------- merge
def merge_partials(out_dir: str, family: str, expected_dirs: int) -> dict:
    out = Path(out_dir)
    parts = [json.loads(p.read_text()) for p in out.glob(f"*_{family}_*.partial")]
    assert parts, f"no partials for family {family} in {out_dir}"
    covered = sum(p["count"] for p in parts)
    assert covered == expected_dirs, f"coverage {covered} != {expected_dirs}"
    merged = ef.merge_block_maxima(parts)
    merged["full_coverage"] = True
    return merged


# ---------------------------------------------------------- subprocess entry point
def _worker_main():
    """argv: node_name family start count ... out_dir  (blocks given as start,count pairs)."""
    node_name = sys.argv[1]
    out_dir = sys.argv[-1]
    fam = sys.argv[2]
    pairs = sys.argv[3:-1]
    blocks = [Block(fam, int(pairs[i]), int(pairs[i + 1]))
              for i in range(0, len(pairs), 2)]
    run_node_jobs(blocks, out_dir, node_name)


if __name__ == "__main__":
    _worker_main()
