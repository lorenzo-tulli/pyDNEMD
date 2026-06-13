"""
Integration tests for DNEMD vector analysis using the example output frames.

These tests load real Cα GRO frames from examples/output/ and check that
collect_vectors + compute_statistics reproduce the known results.

Known values (from examples/output/results/dnemd/summary.csv):
  time_point_ps=1000, n_samples=10, n_ca=313,
  n_ca_significant=266, mean_disp=0.419 Å, max_disp=2.507 Å
"""
import pytest
import numpy as np
from pathlib import Path

from dnemd.dnemd_analysis import collect_vectors, compute_statistics

EXAMPLES = Path(__file__).parent.parent / "examples" / "output"
RUNS = range(1, 3)          # runs 1 and 2
TIME_RANGE_NS = range(1, 6) # 1, 2, 3, 4, 5 ns  →  10 samples total


@pytest.fixture(scope="module")
def vectors_1000ps():
    return collect_vectors(
        base_dir=EXAMPLES,
        time_point_ps=1000,
        runs=RUNS,
        time_range_ns=TIME_RANGE_NS,
    )


@pytest.fixture(scope="module")
def vectors_0ps():
    return collect_vectors(
        base_dir=EXAMPLES,
        time_point_ps=0,
        runs=RUNS,
        time_range_ns=TIME_RANGE_NS,
    )


# ---------------------------------------------------------------------------
# collect_vectors — shape
# ---------------------------------------------------------------------------

def test_collect_vectors_n_samples(vectors_1000ps):
    assert vectors_1000ps.shape[0] == 10   # 2 runs × 5 ns windows


def test_collect_vectors_n_ca(vectors_1000ps):
    assert vectors_1000ps.shape[1] == 313


def test_collect_vectors_xyz(vectors_1000ps):
    assert vectors_1000ps.shape[2] == 3


# ---------------------------------------------------------------------------
# compute_statistics at 1000 ps — regression against known results
# ---------------------------------------------------------------------------

def test_n_ca_significant_1se_1000ps(vectors_1000ps):
    _, _, _, _, adj_disp, _ = compute_statistics(vectors_1000ps, se_threshold=1)
    assert int((adj_disp > 0).sum()) == 266


def test_mean_disp_1000ps(vectors_1000ps):
    _, avg_disp, _, _, _, _ = compute_statistics(vectors_1000ps, se_threshold=1)
    nonzero = avg_disp[avg_disp > 0]
    assert pytest.approx(float(nonzero.mean()), abs=0.005) == 0.419


def test_max_disp_1000ps(vectors_1000ps):
    _, avg_disp, _, _, _, _ = compute_statistics(vectors_1000ps, se_threshold=1)
    assert pytest.approx(float(avg_disp.max()), abs=0.005) == 2.507


# ---------------------------------------------------------------------------
# Sanity check: displacement at t=0 should be much smaller than at t=1000 ps
# (NE and NP start from the same coordinates, diverge only after t=0)
# ---------------------------------------------------------------------------

def test_0ps_displacement_smaller_than_1000ps(vectors_0ps, vectors_1000ps):
    _, avg_disp_0,    _, _, _, _ = compute_statistics(vectors_0ps,    se_threshold=1)
    _, avg_disp_1000, _, _, _, _ = compute_statistics(vectors_1000ps, se_threshold=1)
    assert avg_disp_0.mean() < avg_disp_1000.mean()
