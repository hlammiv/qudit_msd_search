"""Command-line interface: search / reconstruct / asymptotic (NOTES sec 10).

Examples
--------
  python -m qmsd best --p 97 --m 2                # one-liner: best code for a declared (p,m)
  python -m qmsd best --p 101 --m 4 --analytic-only
  python -m qmsd search --p 5 --m 4 --trials 5000
  python -m qmsd search --p 3 --m 4 --sampler capset_climb --target-k 9   # -> [[72,9,3]]_3
  python -m qmsd reconstruct --label "[[20,5,2]]_5"
  python -m qmsd asymptotic --p 5
"""
from __future__ import annotations

import argparse
import sys

# Default trial budget per sampler: the directed samplers do a search *inside* each trial,
# so they need far fewer draws than blind uniform sampling. Keeps the flagship
# `--sampler capset_climb` command fast when the user does not pass --trials.
_DEFAULT_TRIALS = {"uniform": 2000, "capset": 200, "capset_climb": 25, "arc_climb": 15,
                   "plane_spread": 300, "near_cap": 2000, "flat_spread": 300}


def _fmt(c) -> str:
    g = c.gamma
    cert = "" if c.d_certified else "  (d uncertified)"
    ad = "" if c.A_d is None else f"  A_d={c.A_d}"
    return f"{c.label:<22} gamma={g:.4f}{ad}{cert}"


def _cmd_search(args) -> int:
    from .search import search, manhattan_sweep, random_search
    # Resolve the trial budget: an explicit --trials wins; otherwise use the per-sampler
    # default so `--sampler capset_climb` is fast out of the box (its climb makes 2000 draws
    # needlessly slow). The full-scan path (no --m) is uniform, so it gets the uniform default.
    trials = args.trials if args.trials is not None else _DEFAULT_TRIALS[args.sampler]
    if args.m is not None:
        codes = manhattan_sweep(args.p, args.m)
        # Run the randomized search when the space is small, OR whenever the user has asked
        # for a directed search (a non-uniform sampler or a fixed puncture count) -- that is
        # how the paper's d>=3 codes are reached at larger m where uniform sampling stalls.
        directed = args.sampler != "uniform" or args.target_k is not None
        if args.p ** args.m <= 1500 or directed:
            codes += random_search(args.p, args.m, trials, seed=args.seed,
                                   sampler=args.sampler, target_k=args.target_k,
                                   n_jobs=args.jobs, max_triples=args.max_triples,
                                   min_keep_distance=args.min_keep_distance)
        # Preserve structurally distinct puncture sets even when their [[n,k,d]] labels agree;
        # they can have different A_d and geometry and therefore belong in separate artifacts.
        codes = sorted({
            (c.p, c.m, c.n, c.k, c.d, c.A_d, c.w, c.puncture_columns): c
            for c in codes
        }.values(), key=lambda c: c.gamma)
        print(f"# p={args.p}, m={args.m}: {len(codes)} candidate codes (best gamma first)")
        for c in codes[:args.top]:
            print("  " + _fmt(c))
        exported = codes[:args.top]
    else:
        res = search(args.p, trials_per_m=trials, seed=args.seed, top=args.top, n_jobs=args.jobs)
        print(f"# p={args.p}: scanned m={res['scanned']['m_values']}, "
              f"{res['scanned']['n_candidates']} candidates")
        print("## best by yield gamma:")
        for c in res["best_by_gamma"]:
            print("  " + _fmt(c))
        print("## best by single-round cost C (delta_in=1e-3):")
        for c in res["best_by_cost"]:
            print("  " + _fmt(c))
        unique = {}
        for c in res["best_by_gamma"] + res["best_by_cost"]:
            key = (c.p, c.m, c.n, c.k, c.d, c.puncture_columns)
            unique[key] = c
        exported = sorted(unique.values(), key=lambda c: c.gamma)
    if args.output:
        from .catalog import write_search_export
        metadata = {
            "p": args.p, "m": args.m, "trials": trials, "seed": args.seed,
            "top": args.top, "sampler": args.sampler, "target_k": args.target_k,
            "jobs": args.jobs, "max_triples": args.max_triples,
            "min_keep_distance": args.min_keep_distance,
        }
        path = write_search_export(args.output, exported, metadata)
        print(f"# wrote {len(exported)} structured result(s) to {path}")
    return 0


def _cmd_catalog_import(args) -> int:
    from .catalog import import_search_export
    try:
        result = import_search_export(args.source, args.catalog_dir, args.status)
    except (OSError, ValueError, TypeError) as exc:
        print(f"catalog import failed: {exc}", file=sys.stderr)
        return 2
    print(f"catalog: {result['directory']}")
    print(f"imported: {len(result['imported'])}")
    for artifact_id in result["imported"]:
        print(f"  + {artifact_id}")
    print(f"already present: {len(result['unchanged'])}")
    return 0


def _cmd_reconstruct(args) -> int:
    from .oracle import load_oracle
    from .codes import code_from_puncture
    match = [oc for oc in load_oracle() if oc.label.replace(" ", "") == args.label.replace(" ", "")]
    if not match:
        labels = ", ".join(oc.label for oc in load_oracle())
        print(f"no oracle code with label {args.label!r}. Available: {labels}")
        return 1
    from .codes import gamma as _gamma
    oc = match[0]
    # Distance is now certified by the meet-in-the-middle routine even for the large
    # blocks (seconds); A_d exact counting is still only feasible for the small ones.
    c = code_from_puncture(oc.p, oc.m, oc.puncture_columns_1indexed,
                           compute_A_d=(oc.n <= 130), max_distance=oc.d + 1)
    print(f"reconstructed {oc.label} from {oc.k} puncture columns "
          f"(r_max={oc.r_max}, rtilde={oc.r_tilde}):")
    d_show = c.d if c.d is not None else f"{oc.d} (paper value; not certified)"
    cert = "certified" if c.d_certified else "uncertified"
    print(f"  n={c.n}  k={c.k}  d={d_show} ({cert})  full_rank={c.full_rank}  A_d={c.A_d}")
    d_for_gamma = c.d if c.d is not None else oc.d
    print(f"  gamma = log(n/k)/log(d) = {_gamma(c.n, c.k, d_for_gamma):.4f}")
    return 0


def _cmd_asymptotic(args) -> int:
    from .asymptotics import optimal_gamma
    g0, t0 = optimal_gamma(args.p)
    print(f"p={args.p}: asymptotic optimal yield gamma_0 = {g0:.5f} at t_0 = {t0:.5f}")
    return 0


def _cmd_best(args) -> int:
    """One-liner: the best code we can produce for the DECLARED (p, m). Always the closed-form
    Manhattan analytic code; plus the explicit search when p**m is small enough; prints whichever
    has the lower gamma, tagged by provenance."""
    from .search import best_code, EXPLICIT_MAX_BLOCK
    trials = args.trials if args.trials is not None else _DEFAULT_TRIALS.get(args.sampler, 300)
    c = best_code(args.p, args.m, explicit=not args.analytic_only, trials=trials,
                  sampler=args.sampler, target_k=args.target_k, n_jobs=args.jobs, seed=args.seed)
    if c is None:
        print(f"# p={args.p} m={args.m}: no valid code found")
        return 1
    src = ("analytic (Manhattan, closed-form d)" if c.puncture_columns is None
           else "explicit search (exact-certified d<=6)")
    skipped = (not args.analytic_only) and args.p ** args.m > EXPLICIT_MAX_BLOCK
    note = f"  [explicit search skipped: p^m={args.p ** args.m} > {EXPLICIT_MAX_BLOCK}]" if skipped else ""
    print(f"# best code for p={args.p}, m={args.m}  ->  {src}{note}")
    print("  " + _fmt(c))
    return 0


def main(argv=None) -> int:
    """Argparse entry point with search, catalog, reconstruct, asymptotic, and best."""
    parser = argparse.ArgumentParser(
        prog="qmsd",
        description="Discover qudit triorthogonal magic-state-distillation codes "
                    "(after arXiv:2510.10852).",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    ps = sub.add_parser("search", help="search for low-gamma codes given a prime p")
    ps.add_argument("--p", type=int, required=True, help="prime qudit dimension")
    ps.add_argument("--m", type=int, default=None, help="RM variables (default: scan a range)")
    ps.add_argument("--trials", type=int, default=None,
                    help="random-search trials per m (default by sampler: 2000 uniform / "
                         "200 capset / 25 capset_climb -- the climb needs far fewer)")
    ps.add_argument("--seed", type=int, default=0)
    ps.add_argument("--top", type=int, default=10)
    ps.add_argument("--sampler",
                    choices=["uniform", "capset", "capset_climb", "arc_climb", "plane_spread",
                             "near_cap", "flat_spread"],
                    default="uniform",
                    help="puncture sampler (used with --m): 'flat_spread' is the UNIFIED one -- it "
                         "auto-picks the max feasible arc order for the given k (cap / no-4-coplanar "
                         "/ higher), falling back to a near-cap when a cap stalls, and subsumes the "
                         "rest; the fixed operating points are 'capset', 'plane_spread' (no 4 "
                         "coplanar, e.g. [[230,13,6]]_3 d=6), 'near_cap' (+--max-triples, e.g. "
                         "[[200,43,3]]_3); 'arc_climb' climbs the (distance, -A_d) surrogate")
    ps.add_argument("--max-triples", type=int, default=0, dest="max_triples",
                    help="collinear-triple budget for --sampler near_cap (0 = strict cap)")
    ps.add_argument("--target-k", type=int, default=None, dest="target_k",
                    help="fix the puncture count k per candidate (needed for the cap samplers "
                         "to be useful)")
    ps.add_argument("--min-keep-distance", type=int, default=None, dest="min_keep_distance",
                    help="geometric distance FLOOR: skip the (slow) exact-distance MITM on any "
                         "candidate whose fast geometric upper bound (d_RM - max-2flat) is below "
                         "this. Sound; speeds up high-d hunts with cheap samplers (e.g. 6)")
    ps.add_argument("--jobs", type=int, default=1,
                    help="process-level parallelism for the randomized search (trials are "
                         "independent; -1 uses all cores). Serial (1) is bit-for-bit "
                         "reproducible; parallel is reproducible per (seed, sampler, jobs, trials)")
    ps.add_argument("--output", type=str,
                    help="write the displayed candidates as a versioned JSON search bundle")
    ps.set_defaults(func=_cmd_search)

    pc = sub.add_parser("catalog", help="manage explorer catalog artifacts")
    pcsub = pc.add_subparsers(dest="catalog_cmd", required=True)
    pci = pcsub.add_parser("import", help="validate and import a JSON search bundle")
    pci.add_argument("source", help="bundle created by `qmsd search --output ...`")
    pci.add_argument("--catalog-dir", help="override qmsd/data/catalog (also useful in automation)")
    pci.add_argument("--status", choices=["confirmed", "partial", "candidate"],
                     help="override inferred evidence status; confirmed still requires certification")
    pci.set_defaults(func=_cmd_catalog_import)

    pr = sub.add_parser("reconstruct", help="rebuild a published oracle code from its punctures")
    pr.add_argument("--label", type=str, required=True, help='e.g. "[[20,5,2]]_5"')
    pr.set_defaults(func=_cmd_reconstruct)

    pa = sub.add_parser("asymptotic", help="asymptotic optimal yield gamma_0(p)")
    pa.add_argument("--p", type=int, required=True)
    pa.set_defaults(func=_cmd_asymptotic)

    pb = sub.add_parser("best", help="best code for a declared (p,m): analytic Manhattan + explicit "
                                     "search, lower-gamma wins (the one-liner)")
    pb.add_argument("--p", type=int, required=True, help="prime qudit dimension")
    pb.add_argument("--m", type=int, required=True, help="RM variables")
    pb.add_argument("--analytic-only", action="store_true", dest="analytic_only",
                    help="return the closed-form Manhattan code only (skip the explicit search)")
    pb.add_argument("--trials", type=int, default=None,
                    help="explicit-search trials (default by sampler); ignored with --analytic-only")
    pb.add_argument("--sampler",
                    choices=["uniform", "capset", "capset_climb", "arc_climb", "plane_spread",
                             "near_cap", "flat_spread"],
                    default="flat_spread", help="explicit-search sampler (default: flat_spread)")
    pb.add_argument("--target-k", type=int, default=None, dest="target_k",
                    help="fix the puncture count k for the explicit search")
    pb.add_argument("--jobs", type=int, default=1, help="process parallelism for the explicit search")
    pb.add_argument("--seed", type=int, default=0)
    pb.set_defaults(func=_cmd_best)

    args = parser.parse_args(argv)
    return args.func(args)
