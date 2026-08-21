"""Tests for the qmsd CLI (invoked as a subprocess: python -m qmsd ...)."""
import subprocess
import sys
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _run(*args):
    return subprocess.run(
        [sys.executable, "-m", "qmsd", *args],
        cwd=ROOT, capture_output=True, text=True, timeout=120,
    )


def test_cli_reconstruct():
    r = _run("reconstruct", "--label", "[[20,5,2]]_5")
    assert r.returncode == 0, r.stderr
    assert "n=20" in r.stdout and "k=5" in r.stdout and "d=2" in r.stdout
    assert "gamma" in r.stdout


def test_cli_asymptotic():
    r = _run("asymptotic", "--p", "5")
    assert r.returncode == 0, r.stderr
    assert "0.559" in r.stdout  # gamma_0(5) = 0.55914 (ground truth Table 1)


def test_cli_search_smoke():
    r = _run("search", "--p", "5", "--m", "2", "--trials", "40", "--seed", "1")
    assert r.returncode == 0, r.stderr
    assert "gamma=" in r.stdout


def test_cli_search_export_import_and_explorer_discovery(tmp_path, monkeypatch):
    export = tmp_path / "search.json"
    catalog = tmp_path / "catalog"
    searched = _run("search", "--p", "5", "--m", "2", "--trials", "2",
                    "--seed", "4", "--top", "2", "--output", str(export))
    assert searched.returncode == 0, searched.stderr
    bundle = json.loads(export.read_text())
    assert bundle["schema"] == "qmsd.search-results.v1"
    assert bundle["search"]["seed"] == 4
    assert 1 <= len(bundle["codes"]) <= 2
    assert all("distance_certified" in code for code in bundle["codes"])

    imported = _run("catalog", "import", str(export), "--catalog-dir", str(catalog))
    assert imported.returncode == 0, imported.stderr
    assert "imported:" in imported.stdout
    files = sorted(catalog.glob("*.json"))
    assert len(files) == len(bundle["codes"])

    repeated = _run("catalog", "import", str(export), "--catalog-dir", str(catalog))
    assert repeated.returncode == 0, repeated.stderr
    assert f"already present: {len(files)}" in repeated.stdout

    monkeypatch.setenv("QMSD_CATALOG_DIR", str(catalog))
    from qmsd.results import load_result_catalog
    ids = {record.artifact_id for record in load_result_catalog()}
    assert {json.loads(path.read_text())["artifact_id"] for path in files} <= ids


def test_catalog_rejects_false_confirmation(tmp_path):
    bundle = tmp_path / "uncertified.json"
    bundle.write_text(json.dumps({
        "schema": "qmsd.search-results.v1", "search": {}, "codes": [{
            "p": 3, "m": 2, "r": 1, "w": None, "n": 8, "k": 1, "d": None,
            "A_d": None, "gamma": None, "puncture_columns_1indexed": [1],
            "full_rank": True, "distance_certified": False, "A_d_exact": None,
        }]}))
    result = _run("catalog", "import", str(bundle), "--catalog-dir", str(tmp_path / "catalog"),
                  "--status", "confirmed")
    assert result.returncode == 2
    assert "cannot import" in result.stderr


def test_cli_search_capset_climb_reproduces_oracle():
    # The uniform CLI search stalls at [[72,9,2]]; --sampler capset_climb reaches the
    # paper's d=3 code from the command line. No --trials: this exercises the small
    # per-sampler default (25), so the 120s _run timeout also guards against the default
    # regressing back to the slow uniform budget (2000).
    r = _run("search", "--p", "3", "--m", "4", "--sampler", "capset_climb",
             "--target-k", "9", "--seed", "0", "--top", "20")
    assert r.returncode == 0, r.stderr
    assert "[[72,9,3]]_3" in r.stdout
