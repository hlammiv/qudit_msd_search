"""Work partition: split the 3^K message space into independent blocks and assign
them across nodes, emitting per-node job manifests.

Blocks (design sec.6.1): fix the TOP ``t`` message trits -> ``3^t`` independent
blocks, each enumerating the remaining ``K-t`` trits (``3^(K-t)`` messages).  A
block is identified by its integer id in ``[0, 3^t)`` (its base-3 digits are the
fixed values of the block rows).  ``t`` is chosen so a *chunk* of blocks runs in
seconds-to-minutes for fine-grained checkpointing and load balance.

A node's job is a list of CHUNKS; each chunk is a contiguous block range
``(start, count)`` that dwd_core enumerates in one invocation -> one .partial.
The union of all chunks across all nodes is exactly ``[0, 3^t)`` (disjoint,
exhaustive) -- asserted by ``assert_cover``.
"""
from __future__ import annotations

import json
from pathlib import Path


def pow3(e: int) -> int:
    return 3 ** e


def choose_t(K: int, target_blocks: int = 50_000) -> int:
    """Pick t so 3^t >= target_blocks (fine-grained) but t <= K.  ~59049 blocks at K>=10."""
    t = 0
    while pow3(t) < target_blocks and t < K:
        t += 1
    return min(t, K)


def split_contiguous(total_blocks: int, weights: list[float]) -> list[tuple[int, int]]:
    """Partition ``[0, total_blocks)`` into len(weights) contiguous (start, count) ranges
    proportional to ``weights`` (e.g. measured per-node throughput).  Exhaustive, disjoint."""
    assert total_blocks >= 0 and weights and all(w > 0 for w in weights)
    s = sum(weights)
    ranges, start = [], 0
    for i, w in enumerate(weights):
        if i == len(weights) - 1:
            count = total_blocks - start
        else:
            count = int(round(total_blocks * w / s))
            count = max(0, min(count, total_blocks - start))
        ranges.append((start, count))
        start += count
    assert start == total_blocks, (start, total_blocks)
    return ranges


def chunkify(start: int, count: int, chunk_blocks: int) -> list[tuple[int, int]]:
    """Break a (start, count) block range into checkpoint-sized contiguous chunks."""
    assert chunk_blocks >= 1
    out, b = [], start
    end = start + count
    while b < end:
        c = min(chunk_blocks, end - b)
        out.append((b, c))
        b += c
    return out


def make_jobs(K: int, t: int, node_names: list[str], weights: list[float],
              chunk_blocks: int) -> dict:
    """Build the full job plan: per-node list of chunks.  Returns a JSON-able dict."""
    total = pow3(t)
    ranges = split_contiguous(total, weights)
    jobs = {}
    all_chunks = []
    for name, (start, count) in zip(node_names, ranges):
        chunks = chunkify(start, count, chunk_blocks)
        jobs[name] = {
            "block_start": start, "block_count": count,
            "chunks": [{"start": s, "count": c, "chunk_id": f"{name}_{s}_{c}"}
                       for (s, c) in chunks],
        }
        all_chunks += chunks
    assert_cover(all_chunks, total)
    return {
        "K": K, "t": t, "q": 3,
        "total_blocks": total,
        "messages_per_block": pow3(K - t),
        "total_messages": pow3(K),
        "chunk_blocks": chunk_blocks,
        "nodes": jobs,
    }


def assert_cover(chunks: list[tuple[int, int]], total_blocks: int) -> None:
    """Assert the chunks tile [0, total_blocks) exactly: disjoint and exhaustive."""
    covered = sorted(chunks)
    pos = 0
    for (s, c) in covered:
        assert s == pos, f"gap/overlap at block {pos} (chunk starts {s})"
        pos += c
    assert pos == total_blocks, f"coverage {pos} != total {total_blocks}"


def write_plan(plan: dict, path: str | Path) -> None:
    Path(path).write_text(json.dumps(plan, indent=2) + "\n")


def read_plan(path: str | Path) -> dict:
    return json.loads(Path(path).read_text())
