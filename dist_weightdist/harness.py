"""Per-node orchestration: run dwd_core over a node's assigned chunks, verify each
.partial, and record completion in the resume ledger (design sec.6-7).

This is GLUE (not hot): it spawns the compiled core once per chunk.  It is the same
code path a real compute node runs; the two-machine job is just two nodes each
running ``run_node`` over their slice of the plan, writing .partial files into a
shared/collected directory.  No SSH is performed here -- the user copies the binary
+ .g0 to the remote node and runs ``run_node`` there (see RUNBOOK).
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from .checkpoint import Ledger
from .partial import read_partial
from .partition import pow3

DEFAULT_CORE = Path(__file__).resolve().parent / "core_c" / "dwd_core"


def core_binary(explicit: str | None = None) -> Path:
    p = Path(explicit) if explicit else DEFAULT_CORE
    if not p.exists():
        raise FileNotFoundError(
            f"dwd_core binary not found at {p}; build it with core_c/build.sh")
    return p


def run_chunk(core: Path, g0: str, t: int, start: int, count: int,
              out_partial: str, threads: int = 0) -> dict:
    """Invoke dwd_core for one chunk -> .partial; parse + verify it; return its dict."""
    cmd = [str(core), "enum", str(g0), str(t), str(start), str(count), str(out_partial)]
    if threads:
        cmd.append(str(threads))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"dwd_core failed ({res.returncode}): {res.stderr}")
    p = read_partial(out_partial)              # verifies internal checksum
    return p


def run_node(plan: dict, node_name: str, g0: str, out_dir: str,
             ledger_path: str, core: str | None = None, threads: int = 0) -> dict:
    """Run all chunks for ``node_name`` in the plan, resuming from the ledger.

    Each chunk writes ``<out_dir>/<chunk_id>.partial``; completion is recorded only
    after the .partial's internal checksum AND its expected message count are verified.
    Returns a summary {chunks_total, chunks_run, chunks_skipped, messages}.
    """
    K, t = plan["K"], plan["t"]
    mpb = plan["messages_per_block"]
    assert mpb == pow3(K - t)
    coreb = core_binary(core)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ledger = Ledger(ledger_path)
    job = plan["nodes"][node_name]

    run = skipped = 0
    msgs = 0
    for ch in job["chunks"]:
        cid = ch["chunk_id"]
        s, c = ch["start"], ch["count"]
        if c == 0:
            continue
        partial_path = out / f"{cid}.partial"
        expected = c * mpb
        if ledger.is_done(cid) and partial_path.exists():
            skipped += 1
            msgs += expected
            continue
        p = run_chunk(coreb, g0, t, s, c, str(partial_path), threads)
        assert p["checksum"] == expected, (
            f"chunk {cid}: checksum {p['checksum']} != expected {expected}")
        assert p["K"] == K and p["nblocks"] == c
        ledger.mark_done(cid, str(partial_path), p["checksum"], (s, c))
        run += 1
        msgs += expected
    return {"node": node_name, "chunks_total": len(job["chunks"]),
            "chunks_run": run, "chunks_skipped": skipped, "messages": msgs}
