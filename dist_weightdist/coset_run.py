"""Run the ``dwq_coset`` C kernel (A and its cosets A + h_i) and merge into a profile.

Local usage (explicit thread count -- never defaults to all cores):

    python -m dist_weightdist.coset_run <artifact_id> --threads 4 [--outputs 0,5,17] [--t 2]

For a remote box the same pieces split: ``write_cos`` -> scp -> ``dwq_coset enum`` over
disjoint block ranges -> scp back the .cpart files -> ``merge_cparts``.
"""
from __future__ import annotations

import argparse
import struct
import subprocess
from pathlib import Path

import numpy as np

KERNEL = Path(__file__).resolve().parent / "core_c" / "dwq_coset"


def write_cos(path, G0, H1, q):
    G0 = np.asarray(G0, dtype=np.int64) % q
    H1 = np.asarray(H1, dtype=np.int64).reshape(-1, G0.shape[1]) % q
    with open(path, "w") as fh:
        fh.write(f"qmsd-coset v1 {q} {G0.shape[0]} {H1.shape[0]} {G0.shape[1]}\n")
        for r in list(G0) + list(H1):
            fh.write(" ".join(map(str, r.tolist())) + "\n")
    return Path(path)


def read_cpart(path):
    raw = Path(path).read_bytes()
    assert raw[:8] == b"DWQC0001", f"{path}: bad magic"
    q, K, L, n = struct.unpack_from("<4I", raw, 8)
    bs, bc, t = struct.unpack_from("<3Q", raw, 24)
    hist = np.frombuffer(raw, dtype="<i8", offset=48).reshape(1 + L, n + 1)
    return {"q": q, "K": K, "L": L, "n": n, "block_start": bs, "block_count": bc, "t": t,
            "hist": hist}


def merge_cparts(paths):
    """Sum partial histograms; assert the block ranges tile [0, q^t) exactly once."""
    parts = [read_cpart(p) for p in paths]
    q, K, L, n, t = (parts[0][k] for k in ("q", "K", "L", "n", "t"))
    assert all((p["q"], p["K"], p["L"], p["n"], p["t"]) == (q, K, L, n, t) for p in parts)
    ranges = sorted((p["block_start"], p["block_start"] + p["block_count"]) for p in parts)
    pos = 0
    for a, b in ranges:
        assert a == pos, f"block ranges not contiguous at {pos}"
        pos = b
    assert pos == q ** t, "block ranges do not cover all q^t blocks"
    total = sum(p["hist"].astype(object) for p in parts)
    rows = [[int(v) for v in r] for r in total]
    assert all(sum(r) == q ** K for r in rows)
    return {"q": q, "K": K, "A": rows[0], "cosets": rows[1:]}


def run_local(G0, H1, q, workdir, threads, t=2, chunks=1):
    """Enumerate on this machine with ``threads`` OpenMP threads; returns merged histograms."""
    assert threads >= 1
    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    cos = write_cos(work / "in.cos", G0, H1, q)
    nblk = q ** t
    bounds = np.linspace(0, nblk, chunks + 1).astype(int)
    parts = []
    for i in range(chunks):
        out = work / f"part{i}.cpart"
        subprocess.run([str(KERNEL), "enum", str(cos), str(t), str(bounds[i]),
                        str(bounds[i + 1] - bounds[i]), str(out), str(threads)], check=True)
        parts.append(out)
    return merge_cparts(parts)


def main(argv=None):
    from qmsd.performance import profile_from_histograms, record_matrices, summarize
    from qmsd.results import load_result_catalog

    ap = argparse.ArgumentParser()
    ap.add_argument("artifact_id")
    ap.add_argument("--threads", type=int, required=True)
    ap.add_argument("--outputs", default="", help="comma-separated 0-based output indices (default all)")
    ap.add_argument("--t", type=int, default=2)
    ap.add_argument("--workdir", default="runs/perf/coset")
    a = ap.parse_args(argv)
    rec = {r.artifact_id: r for r in load_result_catalog()}[a.artifact_id]
    G0, H1 = record_matrices(rec)
    idx = [int(x) for x in a.outputs.split(",")] if a.outputs else list(range(H1.shape[0]))
    m = run_local(G0, H1[idx], rec.p, Path(a.workdir) / a.artifact_id, a.threads, t=a.t)
    prof = profile_from_histograms(rec.p, m["K"], m["A"], m["cosets"])
    print(rec.label, "outputs", len(idx), summarize(prof, "avg"))


if __name__ == "__main__":
    main()
