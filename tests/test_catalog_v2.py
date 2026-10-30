"""Catalog v2 contract (CATALOG.md): every record carries a construction and a certificate, the
status is what the certificate says, rebuildable records rebuild to their (n, k), stored
witnesses are genuine code words of the stated weight, and v1 files upgrade cleanly."""
import json
import numpy as np
import pytest

from qmsd.catalog import (ACCEPTED_SCHEMAS, CATALOG_SCHEMA, catalog_directory, status_from_certificate,
                          upgrade_record, migrate_catalog)
from qmsd.performance import record_matrices
from qmsd.results import load_result_catalog


def _records():
    return list(load_result_catalog())


def test_every_record_has_construction_and_consistent_status():
    for r in _records():
        assert r.construction is not None and r.construction.get("family") in (
            "punctured_rm", "reed_solomon", "weighted_rm", "manhattan", "none"), r.artifact_id
        if r.status == "candidate":
            assert r.certificate is None or r.certificate.get("lower_bound") is None, r.artifact_id
        else:
            assert r.certificate is not None, r.artifact_id
            assert status_from_certificate(r.certificate) == r.status, r.artifact_id
        if r.certificate and r.certificate.get("exact"):
            assert r.certificate["lower_bound"] == r.certificate["upper_bound"] == r.d, r.artifact_id


def test_on_disk_records_are_v2():
    for path in sorted(catalog_directory().glob("*.json")):
        data = json.loads(path.read_text())
        assert data["schema"] == CATALOG_SCHEMA, path.name
        assert "construction" in data and "certificate" in data and "performance" in data, path.name


def test_rebuild_matches_parameters_and_witnesses_are_code_words():
    """Rebuild every record that fits in memory cheaply; verify shapes and any stored witness."""
    checked = witnesses = 0
    for r in _records():
        fam = r.construction["family"]
        if fam == "none" or (fam == "manhattan" and r.p ** (r.m or 1) > 10 ** 4):
            continue
        if r.p ** (r.m or 1) > 2500 and fam != "weighted_rm":
            continue
        G0, H1 = record_matrices(r)
        assert G0.shape[1] == r.n and H1.shape == (r.k, r.n), r.artifact_id
        checked += 1
        wit = (r.certificate or {}).get("witness")
        if wit:
            vec = np.zeros(r.n, dtype=np.int64)
            vec[wit["columns"]] = wit["coefficients"]
            assert int((vec != 0).sum()) == wit["weight"] == r.certificate["upper_bound"], r.artifact_id
            assert ((G0 @ vec) % r.p == 0).all(), r.artifact_id
            witnesses += 1
    assert checked >= 100 and witnesses >= 5


def test_weighted_records_rebuild_and_certificates_say_what_they_mean():
    recs = {r.artifact_id: r for r in _records()}
    exact = next(r for r in recs.values() if r.label == "[[1756,260,7]]_7")
    assert exact.construction["family"] == "weighted_rm" and exact.status == "confirmed"
    assert exact.certificate["exact"] and exact.certificate["witness"]["weight"] == 7
    G0, H1 = record_matrices(exact)
    assert G0.shape == (66, 1756) and H1.shape == (260, 1756)
    # the two certified k = 270 draws (seeds 19 and 407): both exact d = 7 with stored block/planar witnesses
    k270 = [r for r in recs.values() if r.label == "[[1746,270,7]]_7"]
    assert len(k270) == 2
    for r in k270:
        assert r.status == "confirmed" and r.certificate["exact"] and r.certificate["witness"]["weight"] == 7
        assert r.distance_certified and r.distance_lower_bound == 7
    seed407 = recs["imported-p7-n1746-k270-c833577e4463"]
    assert seed407.geometry["plane_distance"] == 8 and "3+4" in seed407.certificate["upper_method"]


def test_v1_record_upgrades_in_memory_and_on_disk(tmp_path):
    v1 = {"schema": "qmsd.catalog-record.v1", "artifact_id": "x", "label": "[[8,1,2]]_3", "p": 3, "m": 2, "r": 1,
          "n": 8, "k": 1, "d": 2, "puncture_columns": [1], "status": "confirmed", "distance_evidence": "exact", "family": "t",
          "provenance": "t", "source": "t", "note": "", "A_d": None, "gamma_value": None, "w": None, "Ad_evidence": ""}
    up = upgrade_record(v1)
    assert up["schema"] == CATALOG_SCHEMA and up["construction"]["family"] == "punctured_rm"
    assert up["certificate"]["exact"] and status_from_certificate(up["certificate"]) == "confirmed"
    assert upgrade_record(up) == up
    (tmp_path / "x.json").write_text(json.dumps(v1))
    res = migrate_catalog(tmp_path)
    assert res["upgraded"] == ["x.json"]
    assert json.loads((tmp_path / "x.json").read_text())["schema"] == CATALOG_SCHEMA
    assert migrate_catalog(tmp_path)["upgraded"] == []
    assert "qmsd.catalog-record.v1" in ACCEPTED_SCHEMAS
