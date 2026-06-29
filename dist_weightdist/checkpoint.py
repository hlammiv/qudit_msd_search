"""Checkpoint / resume ledger for an hours-long distributed run (design sec.7).

An append-only JSONL ledger records, per completed CHUNK,
``{chunk_id, partial_file, checksum, blocks}``.  A chunk is "done" ONLY after its
.partial is fully written + fsync'd (dwd_core writes .tmp then atomically renames)
and its checksum verified.  On restart the harness reads the ledger, subtracts the
completed chunks, and re-dispatches only the remainder.  Because chunks are
independent and each writes its own checksum'd .partial, resume is exact with no
double counting; a crash mid-chunk leaves no ledger entry, so the chunk is simply
recomputed.

Atomicity: each ledger append is written to a temp line and flushed+fsync'd; the
JSONL format means a torn final line is detected (json parse fails) and dropped on
load, never corrupting earlier entries.
"""
from __future__ import annotations

import json
import os
from pathlib import Path


class Ledger:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.done: dict[str, dict] = {}
        if self.path.exists():
            self._load()

    def _load(self) -> None:
        for line in self.path.read_text().splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                # torn final line from a crash mid-append: ignore (chunk recomputed).
                continue
            self.done[rec["chunk_id"]] = rec

    def is_done(self, chunk_id: str) -> bool:
        return chunk_id in self.done

    def mark_done(self, chunk_id: str, partial_file: str, checksum: int,
                  blocks: tuple[int, int]) -> None:
        rec = {"chunk_id": chunk_id, "partial_file": str(partial_file),
               "checksum": int(checksum), "blocks": list(blocks)}
        with open(self.path, "a") as fp:
            fp.write(json.dumps(rec) + "\n")
            fp.flush()
            os.fsync(fp.fileno())
        self.done[chunk_id] = rec

    def partial_files(self) -> list[str]:
        return [rec["partial_file"] for rec in self.done.values()]

    def total_checksum(self) -> int:
        return sum(int(rec["checksum"]) for rec in self.done.values())
