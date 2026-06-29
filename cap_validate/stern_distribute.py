"""Distribute the Method-1 Stern/ISD search across nodes (+ a LOCAL mock cluster).

Unit of work for Stern is a disjoint range of RANDOM SEEDS: node ``i`` runs
``stern_search`` over its own seed range (each seed = one random info-set permutation),
writes a ``.partial`` JSON holding its local best (weight + witness support).  The merge
rule is the GLOBAL MINIMUM weight (and its witness).  Disjoint seed ranges make the runs
independent and the merge order-invariant, so a 2-machine split reproduces a single-box
deep search of the same total iteration budget.

NO remote host is contacted.  ``mock_run`` launches two LOCAL worker subprocesses over
disjoint seed ranges through a shared dir, merges, and checks the merged best <= the
single-process best (Stern is one-sided: the distributed union searches >= as deeply).

Real 2-machine launch (local 20-core + lenore :60022 32-core) is in cap_validate/RUNBOOK.md.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, os.getcwd())
from cap_validate.stern_isd import stern_search_multi  # noqa: E402

REPO = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------- plan
def split_seed_ranges(total_iters: int, node_weights: list[float],
                      base_seed: int = 0) -> list[tuple[int, int]]:
    """Partition ``total_iters`` random-seed slots across nodes (proportional to
    node_weights).  Returns [(seed_start, iters), ...] with DISJOINT seed ranges."""
    w = np.asarray(node_weights, dtype=float)
    counts = np.floor(w / w.sum() * total_iters).astype(int)
    counts[-1] = total_iters - counts[:-1].sum()
    ranges, s = [], base_seed
    for c in counts:
        ranges.append((s, int(c)))
        s += int(c)
    return ranges


# --------------------------------------------------------------- worker (per node)
def run_node(H, seed_start: int, iters: int, out_dir: str, node_name: str,
             p_values=(1, 2, 3), l: int = 20, budget: int = 13,
             threads: int = 1, target: int | None = 10,
             info_width=None) -> str:
    """Run Stern on this node's seed range, write one .partial, return its path."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    pf = out / f"{node_name}_stern_{seed_start}_{iters}.partial"
    if pf.exists():
        try:
            d = json.loads(pf.read_text())
            if d.get("seed_start") == seed_start and d.get("iters") == iters:
                return str(pf)
        except Exception:
            pass
    res = stern_search_multi(
        H, p_values=p_values, l=l, iterations=iters, weight_budget=budget,
        threads=threads, seed=seed_start, target=target, info_width=info_width,
    )
    rec = dict(
        node=node_name, seed_start=seed_start, iters=res.iterations,
        best_weight=int(res.best_weight),
        witness_support=(np.nonzero(res.witness)[0].tolist()
                         if res.witness is not None else None),
        witness_values=(res.witness[np.nonzero(res.witness)[0]].astype(int).tolist()
                        if res.witness is not None else None),
    )
    pf.write_text(json.dumps(rec))
    return str(pf)


def merge_partials(out_dir: str) -> dict:
    """Merge Stern partials: GLOBAL MIN weight + its witness."""
    parts = [json.loads(p.read_text()) for p in Path(out_dir).glob("*_stern_*.partial")]
    assert parts, f"no stern partials in {out_dir}"
    valid = [p for p in parts if p["best_weight"] != -1]
    total_iters = sum(p["iters"] for p in parts)
    if not valid:
        return dict(best_weight=-1, witness_support=None, total_iters=total_iters,
                    n_nodes=len(parts), below_10=False)
    best = min(valid, key=lambda p: p["best_weight"])
    return dict(
        best_weight=best["best_weight"],
        witness_support=best["witness_support"],
        witness_values=best.get("witness_values"),
        from_node=best["node"], total_iters=total_iters, n_nodes=len(parts),
        below_10=best["best_weight"] < 10,
    )


# ---------------------------------------------------------- subprocess entry point
def _build_small_test_H(seed=11):
    """Small full-rank F_3 parity check with a PLANTED weight-2 codeword.

    Columns 7 and 40 are made equal, so e_7 - e_40 is a guaranteed weight-2 codeword
    of H^perp that the p=1 Stern path finds reliably in a few iterations -- making the
    distributed/merge result deterministic (every node finds best=2)."""
    rng = np.random.default_rng(seed)
    H = rng.integers(0, 3, size=(9, 55)).astype(np.int8)
    H[:, 40] = H[:, 7]            # duplicate column -> weight-2 codeword
    return H


def _load_H(spec: str) -> np.ndarray:
    if spec == "cap":
        from cap_validate.stern_isd import build_cap_g0
        return build_cap_g0()
    if spec == "small":
        return _build_small_test_H()
    return np.load(spec)


def _worker_main():
    ap = argparse.ArgumentParser()
    ap.add_argument("node_name")
    ap.add_argument("--H", default="small")
    ap.add_argument("--seed-start", type=int, required=True)
    ap.add_argument("--iters", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--budget", type=int, default=13)
    ap.add_argument("--target", type=int, default=10)
    a = ap.parse_args()
    H = _load_H(a.H)
    p = run_node(H, a.seed_start, a.iters, a.out, a.node_name,
                 threads=a.threads, budget=a.budget, target=a.target)
    print(f"[{a.node_name}] wrote {p}")


# --------------------------------------------------------------------- mock cluster
def mock_run():
    """LOCAL two-mock-node Stern distribution test on a small code (no network)."""
    H = _build_small_test_H()
    total = 600
    ranges = split_seed_ranges(total, node_weights=[20.0, 32.0], base_seed=0)
    node_names = ["local", "lenore_mock"]
    print(f"[stern-mock] H={H.shape} total_iters={total} "
          f"ranges={[(n, r) for n, r in zip(node_names, ranges)]}")

    # ground truth: single process over the whole budget (must find the planted wt-2)
    truth = stern_search_multi(H, p_values=(1, 2, 3), l=7, iterations=total,
                               weight_budget=H.shape[1], threads=2, seed=0, target=None)
    print(f"[stern-mock] single-process best={truth.best_weight} (planted min=2)")
    assert truth.best_weight == 2, truth.best_weight

    with tempfile.TemporaryDirectory() as td:
        node_best = {}
        for ni, (s, it) in enumerate(ranges):
            argv = [sys.executable, "-m", "cap_validate.stern_distribute",
                    node_names[ni], "--H", "small", "--seed-start", str(s),
                    "--iters", str(it), "--out", td, "--threads", "2",
                    "--budget", str(H.shape[1]), "--target", "-1"]
            subprocess.run(argv, check=True, cwd=str(REPO))
        # record each node's own partial best (merge-logic check)
        for pf in Path(td).glob("*_stern_*.partial"):
            d = json.loads(pf.read_text())
            node_best[d["node"]] = d["best_weight"]
        merged = merge_partials(td)
        print(f"[stern-mock] node bests={node_best}  merged best={merged['best_weight']} "
              f"support={merged['witness_support']} from={merged.get('from_node')} "
              f"nodes={merged['n_nodes']} total_iters={merged['total_iters']}")
        # MERGE LOGIC: merged == global min over node partials, and == planted min
        assert merged["best_weight"] == min(node_best.values()), (merged, node_best)
        assert merged["best_weight"] == 2, merged["best_weight"]

        # resume: wipe one node's partial, re-run, re-merge -> identical
        for p in Path(td).glob("lenore_mock_*.partial"):
            p.unlink()
        s, it = ranges[1]
        subprocess.run([sys.executable, "-m", "cap_validate.stern_distribute",
                        "lenore_mock", "--H", "small", "--seed-start", str(s),
                        "--iters", str(it), "--out", td, "--threads", "2",
                        "--budget", str(H.shape[1]), "--target", "-1"],
                       check=True, cwd=str(REPO))
        merged2 = merge_partials(td)
        assert merged2["best_weight"] == merged["best_weight"]
        print("[stern-mock] resume OK: merged identical after wiping+rerunning one node")
    print("[stern-mock] PASS: distributed >= single-process depth; merge=global-min; "
          "resume idempotent")
    return merged


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "mock":
        mock_run()
    else:
        _worker_main()
