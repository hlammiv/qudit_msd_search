"""Exact distillation performance engine vs arXiv:2610.08012 Table 2 and the repo's A_d."""
import pytest

from qmsd.cyclic import construction_A, construction_B
from qmsd.mc_distill import simulate
from qmsd.performance import (canonical_logical_rows, leading_order, profile_exact,
                              puncture_columns_of, record_matrices, threshold)
from qmsd.reedmuller import r_max, rm_generator
from qmsd.results import load_result_catalog
from qmsd.triorthogonal import build_triorthogonal_code


def _rm(p, m):
    r = r_max(m, p)
    G = rm_generator(r, m, p)
    G0 = build_triorthogonal_code(p, m, r, [1], G=G)["G0"]
    return G0, canonical_logical_rows(G, [0], p)


# (name, builder, p, d, c_eps, threshold, P_acc(0.01))  -- paper Table 2
TABLE2 = [
    ("[[19,1,3]]3 B", lambda: construction_B(3, 20, {2, 5, 6, 14, 15, 18}), 3, 3, 5.25, 0.226, 0.83),
    ("[[40,1,4]]3 A", lambda: construction_A(3, 40, {1, 2, 3, 5, 6, 9, 14, 15, 18, 27}), 3, 4, 16.25, 0.223, 0.67),
    ("[[8,1,2]]3 RM", lambda: _rm(3, 2), 3, 2, 2.0, 0.211, 0.92),
    ("[[26,1,2]]3 RM", lambda: _rm(3, 3), 3, 2, 6.5, 0.066, 0.77),
    ("[[11,1,3]]5 B", lambda: construction_B(5, 12, {2, 3, 9, 10}), 5, 3, 1.44, 0.346, 0.90),
    ("[[24,1,3]]5 RM", lambda: _rm(5, 2), 5, 3, 1.5, 0.263, 0.79),
    ("[[4,1,2]]5 RS", lambda: _rm(5, 1), 5, 2, 1.5, 0.363, 0.96),
    ("[[16,1,4]]5", lambda: construction_A(5, 16, {1, 4, 5, 9, 13}), 5, 4, 3.38, 0.363, 0.85),
    ("[[7,1,3]]7 B", lambda: construction_B(7, 8, {1, 4, 7}), 7, 3, 0.97, 0.446, 0.93),
    ("[[12,1,3]]7 A", lambda: construction_A(7, 12, {1, 2, 7}), 7, 3, 1.11, 0.366, 0.89),
    ("[[6,1,2]]7 RS", lambda: _rm(7, 1), 7, 2, 2.5, 0.232, 0.94),
]


@pytest.mark.parametrize("name,build,p,d,c,thr,pacc", TABLE2, ids=[t[0] for t in TABLE2])
def test_paper_table2(name, build, p, d, c, thr, pacc):
    prof = profile_exact(*build(), p)
    lo = leading_order(prof, "block")
    assert lo.d == d
    assert lo.c_eps == pytest.approx(c, abs=0.006)
    assert threshold(prof, "block").eps == pytest.approx(thr, abs=0.0006)
    assert float(prof.p_accept(0.01)) == pytest.approx(pacc, abs=0.006)
    # k = 1: per-output and block coincide (same function, different integer normalisation)
    for eps in (1e-3, 0.05, 0.2):
        assert float(prof.eps_out(eps, "avg")) == pytest.approx(float(prof.eps_out(eps, "block")), rel=1e-12)


def _record(aid):
    return {r.artifact_id: r for r in load_result_catalog()}[aid]


@pytest.mark.parametrize("aid", ["paper-table3-p5-m2-n20-k5", "paper-table3-p3-m4-n72-k9",
                                 "paper-table3-p5-m2-n24-k1"])
def test_block_leading_order_reproduces_catalog_Ad(aid):
    rec = _record(aid)
    prof = profile_exact(*record_matrices(rec), rec.p)
    lo = leading_order(prof, "block")
    assert (lo.d, lo.L_d) == (rec.d, rec.A_d)
    # per-output count sums over outputs, so each output sees at least A_d/k
    assert leading_order(prof, "avg").L_d >= lo.L_d / rec.k


def test_small_manhattan_records_rebuild_explicitly():
    # the [[17,6,3]]_23 Table-2 code is small enough to rebuild from its weight-<=w point set
    rec = _record("manhattan-table2-p23-m1-w5")
    assert len(puncture_columns_of(rec)) == rec.k
    lo = leading_order(profile_exact(*record_matrices(rec), rec.p, kinds=("avg",)), "avg")
    assert lo.d == rec.d
    # astronomically large Manhattan codes stay analytic-only
    assert puncture_columns_of(_record("manhattan-table2-p5-m16-w16")) is None


def test_monte_carlo_agrees():
    rec = _record("paper-table3-p5-m2-n20-k5")
    G0, H1 = record_matrices(rec)
    prof = profile_exact(G0, H1, 5)
    mc = simulate(G0, H1, 5, 0.03, shots=1_000_000, seed=7)
    assert abs(mc["P_acc"] - float(prof.p_accept(0.03))) < 4 * mc["P_acc_err"]
    assert abs(mc["eps_out_block"] - float(prof.eps_out(0.03, "block"))) < 4 * mc["eps_out_block_err"]
    assert abs(mc["eps_out_avg"] - float(prof.eps_out(0.03, "avg"))) < 4 * mc["eps_out_avg_err"]


def test_dwq_coset_kernel_matches_numpy(tmp_path):
    from dist_weightdist.coset_run import KERNEL, run_local
    from qmsd.performance import profile_from_histograms
    from qmsd.weightdist import coset_weight_enumerators

    if not KERNEL.exists():
        pytest.skip("dwq_coset not built (gcc -O3 -march=native -fopenmp dwq_coset.c)")
    rec = _record("optimized-p5-m3-ad396")
    G0, H1 = record_matrices(rec)
    m = run_local(G0, H1, rec.p, tmp_path, threads=1, t=2, chunks=2)
    py = coset_weight_enumerators(G0, [None] + list(H1), rec.p)
    assert m["A"] == py[0] and m["cosets"] == py[1:]
    # output subsets average the per-output polynomials over the chosen outputs only
    sub = profile_from_histograms(rec.p, m["K"], m["A"], m["cosets"][:3])
    assert sub.k == 3
