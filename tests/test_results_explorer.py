"""Tests for the normalized research catalog and Streamlit explorer."""
import math
from pathlib import Path

import numpy as np
import pytest

from qmsd.results import (
    affine_line_profile, distillation_series, code_structure, load_result_catalog,
    pareto_front, puncture_points, record_dict,
)


def by_id(artifact_id):
    return next(r for r in load_result_catalog() if r.artifact_id == artifact_id)


def test_catalog_is_curated_and_ids_are_stable():
    records = load_result_catalog()
    assert len(records) >= 31  # 31 built-ins plus any locally imported research artifacts
    assert len({r.artifact_id for r in records}) == len(records)
    assert {r.status for r in records} >= {"confirmed", "candidate", "refuted"}
    assert all(r.n > 0 and r.k > 0 and r.p >= 2 for r in records)


def test_same_code_parameters_preserve_distinct_artifacts():
    paper = by_id("paper-table3-p5-m4-n519-k106")
    optimized = by_id("optimized-p5-flagship-ad1904")
    assert paper.label == optimized.label == "[[519,106,5]]_5"
    assert paper.artifact_id != optimized.artifact_id
    assert paper.A_d == 2180
    assert optimized.A_d == 1904
    assert paper.puncture_columns != optimized.puncture_columns


def test_evidence_policy_is_explicit():
    refuted = by_id("refuted-qutrit-cap")
    candidate = by_id("candidate-p17-d7")
    assert refuted.status == "refuted"
    assert candidate.status == "candidate"
    assert refuted.distance_certified is False
    assert candidate.d is None


def test_distillation_series_requires_confirmed_evidence():
    code = by_id("optimized-qutrit-ad572")
    rows = distillation_series(code, [1e-2, 1e-3])
    assert rows[0]["cost"] > 0
    assert rows[0]["delta_out"] > rows[1]["delta_out"] > 0
    refuted = distillation_series(by_id("refuted-qutrit-cap"), [1e-2])
    assert refuted[0]["delta_out"] is None


def test_distillation_series_handles_astronomical_blocks():
    code = by_id("manhattan-table2-p2-m58-w14")
    row = distillation_series(code, [1e-2])[0]
    assert row["accepted_outputs"] == 0
    assert math.isinf(row["cost"])
    assert row["delta_out"] is None


def test_pareto_excludes_uncertain_claims():
    records = load_result_catalog()
    front = pareto_front(records, "gamma", "n")
    assert front
    assert all(by_id(i).status == "confirmed" for i in front)
    assert "refuted-qutrit-cap" not in front


def test_puncture_geometry_and_serialization():
    code = by_id("paper-table3-p3-m4-n72-k9")
    points = puncture_points(code)
    assert len(points) == code.k
    assert all(len(point) == code.m for point in points)
    data = record_dict(code)
    assert data["gamma"] == pytest.approx(code.gamma)
    assert data["rate"] == pytest.approx(code.k / code.n)


def test_code_structure_maps_logical_to_physical_qudits():
    code = by_id("paper-table3-p3-m4-n72-k9")
    structure = code_structure(code)
    assert structure is not None and structure["full_rank"]
    assert structure["X_stab"].shape[1] == code.n
    assert structure["Z_stab"].shape[1] == code.n
    assert structure["logical_X"].shape == (code.k, code.n)
    assert len(structure["physical_columns"]) == code.n
    # Logical-X representatives commute with both CSS stabilizer groups.
    assert np.all(structure["X_stab"] @ structure["logical_X"].T % code.p == 0)
    assert np.all(structure["Z_stab"] @ structure["logical_X"].T % code.p == 0)


def test_affine_line_profile_partitions_small_geometry():
    code = by_id("paper-table3-p3-m4-n72-k9")
    profile = affine_line_profile(code)
    assert profile is not None
    expected = code.p ** (code.m - 1) * (code.p ** code.m - 1) // (code.p - 1)
    assert profile["line_count"] == expected
    assert sum(count for _, count in profile["histogram"]) == expected
    assert len(profile["max_line_columns"]) == code.p
    assert profile["max_occupancy"] <= code.p


def test_streamlit_app_smoke():
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    app_path = Path(__file__).resolve().parents[1] / "app.py"
    app = AppTest.from_file(app_path, default_timeout=30).run()
    assert not app.exception
    assert app.title[0].value == "QMSD Research Explorer"
    assert len(app.dataframe) >= 1
    assert any(x.label == "Inspect artifact" for x in app.selectbox)
