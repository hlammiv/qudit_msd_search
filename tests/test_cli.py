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


def test_cli_search_capset_climb_reproduces_oracle():
    # The uniform CLI search stalls at [[72,9,2]]; --sampler capset_climb reaches the
    # paper's d=3 code from the command line. No --trials: this exercises the small
    # per-sampler default (25), so the 120s _run timeout also guards against the default
    # regressing back to the slow uniform budget (2000).
    r = _run("search", "--p", "3", "--m", "4", "--sampler", "capset_climb",
             "--target-k", "9", "--seed", "0", "--top", "20")
    assert r.returncode == 0, r.stderr
    assert "[[72,9,3]]_3" in r.stdout
