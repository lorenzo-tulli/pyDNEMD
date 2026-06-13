"""
Integration tests for equilibrium analysis using the example XVG files.

These tests load real RMSD/RMSF XVG files from examples/output/results/equilibrium/
and check that parse_xvg + the summary logic reproduce the known results.

Known values (from examples/output/results/equilibrium/summary.csv):
  run 1: mean_rmsd=1.628 Å, max_rmsd=2.082 Å, mean_rmsf=0.847 Å, max_rmsf=4.787 Å
  run 2: mean_rmsd=1.407 Å, max_rmsd=1.815 Å, mean_rmsf=0.766 Å, max_rmsf=3.268 Å
"""
import pytest
import numpy as np
from pathlib import Path

from dnemd.analysis import parse_xvg

RESULTS = Path(__file__).parent.parent / "examples" / "output" / "results" / "equilibrium"


def summarise_xvg(xvg_path: Path) -> dict:
    """Replicate EquilibriumAnalyser._summarise for one XVG (RMSD or RMSF)."""
    _, values = parse_xvg(xvg_path)
    values_A = values * 10   # nm -> Å
    return {
        "mean": round(float(np.mean(values_A)), 3),
        "max":  round(float(np.max(values_A)),  3),
    }


# ---------------------------------------------------------------------------
# parse_xvg — basic sanity
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("filename", [
    "rmsd_run1.xvg", "rmsd_run2.xvg",
    "rmsf_run1.xvg", "rmsf_run2.xvg",
])
def test_parse_xvg_nonempty(filename):
    x, y = parse_xvg(RESULTS / filename)
    assert len(x) > 0
    assert len(x) == len(y)


@pytest.mark.parametrize("filename", [
    "rmsd_run1.xvg", "rmsd_run2.xvg",
    "rmsf_run1.xvg", "rmsf_run2.xvg",
])
def test_parse_xvg_positive_values(filename):
    _, y = parse_xvg(RESULTS / filename)
    assert (y > 0).all()


# ---------------------------------------------------------------------------
# RMSD summary — regression against known results
# ---------------------------------------------------------------------------

def test_rmsd_mean_run1():
    s = summarise_xvg(RESULTS / "rmsd_run1.xvg")
    assert pytest.approx(s["mean"], abs=0.005) == 1.628


def test_rmsd_max_run1():
    s = summarise_xvg(RESULTS / "rmsd_run1.xvg")
    assert pytest.approx(s["max"], abs=0.005) == 2.082


def test_rmsd_mean_run2():
    s = summarise_xvg(RESULTS / "rmsd_run2.xvg")
    assert pytest.approx(s["mean"], abs=0.005) == 1.407


def test_rmsd_max_run2():
    s = summarise_xvg(RESULTS / "rmsd_run2.xvg")
    assert pytest.approx(s["max"], abs=0.005) == 1.815


# ---------------------------------------------------------------------------
# RMSF summary — regression against known results
# ---------------------------------------------------------------------------

def test_rmsf_mean_run1():
    s = summarise_xvg(RESULTS / "rmsf_run1.xvg")
    assert pytest.approx(s["mean"], abs=0.005) == 0.847


def test_rmsf_max_run1():
    s = summarise_xvg(RESULTS / "rmsf_run1.xvg")
    assert pytest.approx(s["max"], abs=0.005) == 4.787


def test_rmsf_mean_run2():
    s = summarise_xvg(RESULTS / "rmsf_run2.xvg")
    assert pytest.approx(s["mean"], abs=0.005) == 0.766


def test_rmsf_max_run2():
    s = summarise_xvg(RESULTS / "rmsf_run2.xvg")
    assert pytest.approx(s["max"], abs=0.005) == 3.268
