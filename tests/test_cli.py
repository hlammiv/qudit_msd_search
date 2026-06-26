"""Tests for the qmsd CLI (invoked as a subprocess: python -m qmsd ...)."""
import subprocess
import sys
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
