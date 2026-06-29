"""Reader/locator for the language-agnostic .partial format written by dwd_core.

Byte layout (little-endian) -- frozen in DISTRIBUTED_WEIGHTDIST_DESIGN.md sec.3.2:

    magic   : 8 bytes  b"DWDP0001"
    q       : u32   (=3)
    K       : u32
    n       : u32
    nblocks : u32      number of message-blocks folded into this partial
    checksum: u64      sum_w hist[w]  == total messages covered
    hist    : (n+1) x i64   the int64 weight enumerator for this partial
    manifest: nblocks x u32  the block ids included

The harness verifies ``checksum == sum(hist) == sum(block sizes)`` before accepting.
"""
from __future__ import annotations

import struct
from pathlib import Path

MAGIC = b"DWDP0001"


def read_partial(path: str | Path) -> dict:
    """Parse a .partial file; verify magic + internal checksum (== sum(hist))."""
    data = Path(path).read_bytes()
    assert data[:8] == MAGIC, f"{path}: bad magic {data[:8]!r}"
    q, K, n, nblocks = struct.unpack_from("<IIII", data, 8)
    (checksum,) = struct.unpack_from("<Q", data, 24)
    off = 32
    hist = list(struct.unpack_from(f"<{n + 1}q", data, off))
    off += (n + 1) * 8
    blocks = list(struct.unpack_from(f"<{nblocks}I", data, off))
    off += nblocks * 4
    assert off == len(data), f"{path}: trailing bytes ({len(data) - off})"
    s = sum(hist)
    assert s == checksum, f"{path}: internal checksum {checksum} != sum(hist) {s}"
    assert all(h >= 0 for h in hist), f"{path}: negative histogram entry"
    assert len(blocks) == nblocks
    return {"q": q, "K": K, "n": n, "nblocks": nblocks,
            "checksum": checksum, "hist": hist, "blocks": blocks}
