"""Serialization of the generator G0 for the distributed engine.

The hot enumeration core may be C / Rust / numba; the file format is therefore deliberately
trivial to parse from any language and carries no framework dependency.

FORMAT  (".g0", plain ASCII, LF-terminated lines)
------------------------------------------------------------------------------------------
    Line 1:   "qmsd-g0 v1 <q> <K> <n>"            (magic, version, field, #rows, #cols)
    Lines 2..K+1: each a string of exactly n characters in {0,1,...,q-1}  (one generator row).

So row i, column j is the j-th character of line i+1, interpreted as a base-q digit.  No
spaces inside a row -> a length-n codeword/row is one token; the compiled core can mmap and
index directly.  ``G0`` is K x n over F_q; the engine enumerates all q^K messages m and
histograms wt(m @ G0).

This module reads/writes that format and reduces entries mod q.  ``A_w <= q^K`` (the int64
overflow proof) only needs K and q, both in the header.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

MAGIC = "qmsd-g0"
VERSION = "v1"


def write_g0(path: str | Path, G0, q: int) -> None:
    """Write a K x n F_q generator matrix to ``path`` in the .g0 format."""
    M = np.asarray(G0)
    assert M.ndim == 2, "G0 must be 2D"
    M = (M.astype(np.int64) % q).astype(np.int64)
    K, n = M.shape
    assert q <= 10, "single-digit format supports q <= 10; widen for larger q"
    lines = [f"{MAGIC} {VERSION} {q} {K} {n}"]
    for i in range(K):
        lines.append("".join(str(int(v)) for v in M[i]))
    Path(path).write_text("\n".join(lines) + "\n")


def read_g0(path: str | Path):
    """Read a .g0 file -> (G0 int64 array K x n reduced mod q, q).  Validates the header/shape."""
    raw = Path(path).read_text().splitlines()
    assert raw, "empty .g0 file"
    hdr = raw[0].split()
    assert len(hdr) == 5 and hdr[0] == MAGIC and hdr[1] == VERSION, f"bad header: {raw[0]!r}"
    q, K, n = int(hdr[2]), int(hdr[3]), int(hdr[4])
    body = raw[1 : 1 + K]
    assert len(body) == K, f"expected {K} rows, found {len(body)}"
    M = np.empty((K, n), dtype=np.int64)
    for i, line in enumerate(body):
        assert len(line) == n, f"row {i} has length {len(line)}, expected n={n}"
        row = np.frombuffer(line.encode("ascii"), dtype=np.uint8).astype(np.int64) - ord("0")
        assert np.all((row >= 0) & (row < q)), f"row {i} has a digit outside [0,{q})"
        M[i] = row
    return M % q, q
