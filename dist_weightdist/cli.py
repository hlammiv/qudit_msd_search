"""Command-line entry points for the distributed weight-distribution engine.

  python -m dist_weightdist.cli run    --g0 FILE --out DIR [--threads N] [--t T]
        Single-machine end to end: enumerate all 3^K codewords of G0, merge, run
        the bignum MacWilliams, write A.txt / B.txt / result.json (d, A_d, B).

  python -m dist_weightdist.cli plan   --g0 FILE --out plan.json
                                       --nodes a,b --weights 1,1 [--t T] [--chunk C]
        Emit a 2-node (or N-node) job plan partitioning the 3^t blocks.

  python -m dist_weightdist.cli node   --plan plan.json --node NAME --g0 FILE
                                       --out DIR --ledger L [--threads N]
        Run one node's chunks (resumable via the ledger).  RUN THIS ON EACH MACHINE.

  python -m dist_weightdist.cli merge  --partials DIR/glob ... --g0 FILE --out DIR
        Collect .partial files from all nodes -> certify (d, A_d, B).
"""
from __future__ import annotations

import argparse
import glob
import json
import time
from pathlib import Path

from . import correctness as cx
from . import merge as mg
from . import partition as pt
from .harness import core_binary, run_chunk, run_node
from .io_g0 import read_g0


def _K_of(g0: str) -> tuple[int, int]:
    G, q = read_g0(g0)
    assert q == 3, "engine is q=3 only"
    return G.shape[0], G.shape[1]


def cmd_run(args):
    K, n = _K_of(args.g0)
    cx.assert_int64_safe(K, 3)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    # pick t so there are plenty of blocks to feed the thread pool (>= 8x threads)
    threads = args.threads or 0
    if args.t is not None:
        t = args.t
    else:
        target = max(8 * (threads or 16), 64)
        t = min(pt.choose_t(K, target), K)
    core = core_binary(args.core)
    partial = out / "single.partial"
    t0 = time.time()
    p = run_chunk(core, args.g0, t, 0, pt.pow3(t), str(partial), threads)
    dt = time.time() - t0
    assert p["checksum"] == 3 ** K, f"checksum {p['checksum']} != 3^{K}"
    res = mg.certify([str(partial)], search_kmax=args.search_kmax)
    summary = mg.write_outputs(res, str(out))
    rate = (3 ** K) / dt if dt > 0 else 0
    print(json.dumps({**summary, "enum_seconds": round(dt, 2),
                      "codewords_per_s": round(rate)}, indent=2))
    print(f"d={res['d']}  A_d={res['A_d']}  (n={n}, K={K})  "
          f"{3**K:.3e} codewords in {dt:.1f}s = {rate:.3e} cw/s")


def cmd_plan(args):
    K, n = _K_of(args.g0)
    cx.assert_int64_safe(K, 3)
    nodes = args.nodes.split(",")
    weights = [float(w) for w in args.weights.split(",")] if args.weights else [1.0] * len(nodes)
    assert len(weights) == len(nodes)
    t = args.t if args.t is not None else pt.choose_t(K, args.target_blocks)
    plan = pt.make_jobs(K, t, nodes, weights, args.chunk)
    pt.write_plan(plan, args.out)
    print(json.dumps({"out": args.out, "K": K, "n": n, "t": t,
                      "total_blocks": plan["total_blocks"],
                      "messages_per_block": plan["messages_per_block"],
                      "nodes": {nm: (jb["block_start"], jb["block_count"],
                                     len(jb["chunks"]))
                                for nm, jb in plan["nodes"].items()}}, indent=2))


def cmd_node(args):
    plan = pt.read_plan(args.plan)
    summ = run_node(plan, args.node, args.g0, args.out, args.ledger,
                    core=args.core, threads=args.threads or 0)
    print(json.dumps(summ, indent=2))


def cmd_merge(args):
    paths = []
    for pat in args.partials:
        paths += sorted(glob.glob(pat)) if any(c in pat for c in "*?[") else [pat]
    paths = [p for p in paths if p.endswith(".partial")] or paths
    assert paths, "no .partial files matched"
    res = mg.certify(paths, search_kmax=args.search_kmax)
    summary = mg.write_outputs(res, args.out)
    print(json.dumps({**summary, "n_partials": len(paths)}, indent=2))
    print(f"d={res['d']}  A_d={res['A_d']}")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="dist_weightdist.cli")
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="single-machine end to end")
    r.add_argument("--g0", required=True); r.add_argument("--out", required=True)
    r.add_argument("--threads", type=int, default=0); r.add_argument("--t", type=int, default=None)
    r.add_argument("--core", default=None); r.add_argument("--search-kmax", type=int, default=None)
    r.set_defaults(fn=cmd_run)

    p = sub.add_parser("plan", help="emit an N-node job plan")
    p.add_argument("--g0", required=True); p.add_argument("--out", required=True)
    p.add_argument("--nodes", required=True); p.add_argument("--weights", default=None)
    p.add_argument("--t", type=int, default=None); p.add_argument("--chunk", type=int, default=64)
    p.add_argument("--target-blocks", type=int, default=50_000)
    p.set_defaults(fn=cmd_plan)

    nd = sub.add_parser("node", help="run one node's chunks (resumable)")
    nd.add_argument("--plan", required=True); nd.add_argument("--node", required=True)
    nd.add_argument("--g0", required=True); nd.add_argument("--out", required=True)
    nd.add_argument("--ledger", required=True); nd.add_argument("--threads", type=int, default=0)
    nd.add_argument("--core", default=None)
    nd.set_defaults(fn=cmd_node)

    m = sub.add_parser("merge", help="collect partials -> certify (d, A_d, B)")
    m.add_argument("--partials", nargs="+", required=True); m.add_argument("--out", required=True)
    m.add_argument("--g0", default=None); m.add_argument("--search-kmax", type=int, default=None)
    m.set_defaults(fn=cmd_merge)

    args = ap.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
