import numpy as np
import pytest
from dnemd.dnemd_analysis import compute_statistics, write_stats_txt


def make_vectors(n_samples=20, n_ca=10, seed=42):
    return np.random.default_rng(seed).standard_normal((n_samples, n_ca, 3))


# ---------------------------------------------------------------------------
# compute_statistics — output shapes
# ---------------------------------------------------------------------------

def test_output_shapes():
    vectors = make_vectors(n_samples=20, n_ca=10)
    avg_vectors, avg_disp, se, se_mag, adj_disp, adj_vectors = compute_statistics(vectors)
    assert avg_vectors.shape == (10, 3)
    assert avg_disp.shape  == (10,)
    assert se.shape        == (10, 3)
    assert se_mag.shape    == (10,)
    assert adj_disp.shape  == (10,)
    assert adj_vectors.shape == (10, 3)


def test_adjusted_disp_nonnegative():
    vectors = make_vectors()
    _, _, _, _, adj_disp, _ = compute_statistics(vectors)
    assert (adj_disp >= 0).all()


def test_adjusted_vectors_no_larger_than_original():
    vectors = make_vectors()
    _, avg_disp, _, _, adj_disp, _ = compute_statistics(vectors)
    assert (adj_disp <= avg_disp + 1e-12).all()


# ---------------------------------------------------------------------------
# compute_statistics — se_threshold behaviour
# ---------------------------------------------------------------------------

def test_se_threshold_2_keeps_fewer_atoms():
    # Strong consistent signal so some atoms survive both thresholds
    rng = np.random.default_rng(0)
    vectors = np.ones((50, 20, 3)) + rng.standard_normal((50, 20, 3)) * 0.05
    _, _, _, _, adj_1, _ = compute_statistics(vectors, se_threshold=1)
    _, _, _, _, adj_2, _ = compute_statistics(vectors, se_threshold=2)
    assert (adj_2 > 0).sum() <= (adj_1 > 0).sum()


def test_se_threshold_default_is_1():
    vectors = make_vectors()
    result_default = compute_statistics(vectors)
    result_explicit = compute_statistics(vectors, se_threshold=1)
    np.testing.assert_array_equal(result_default[4], result_explicit[4])


def test_zero_mean_signal_all_filtered():
    # Perfectly symmetric: mean = 0, so every atom should be zeroed out
    half = np.ones((10, 5, 3))
    vectors = np.concatenate([half, -half], axis=0)
    _, _, _, _, adj_disp, _ = compute_statistics(vectors)
    np.testing.assert_array_equal(adj_disp, np.zeros(5))


# ---------------------------------------------------------------------------
# write_stats_txt
# ---------------------------------------------------------------------------

def test_write_stats_txt_creates_file(tmp_path):
    n_ca = 5
    rng = np.random.default_rng(1)
    avg_vectors = rng.standard_normal((n_ca, 3))
    avg_disp    = np.linalg.norm(avg_vectors, axis=1)
    se          = np.abs(rng.standard_normal((n_ca, 3))) * 0.01
    se_mag      = np.abs(rng.standard_normal(n_ca)) * 0.01
    out = tmp_path / "stats.txt"
    write_stats_txt(out, avg_vectors, avg_disp, se, se_mag, n_samples=20)
    assert out.exists()


def test_write_stats_txt_line_count(tmp_path):
    n_ca = 7
    rng = np.random.default_rng(2)
    avg_vectors = rng.standard_normal((n_ca, 3))
    avg_disp    = np.linalg.norm(avg_vectors, axis=1)
    se          = np.ones((n_ca, 3)) * 0.01
    se_mag      = np.ones(n_ca) * 0.01
    out = tmp_path / "stats.txt"
    write_stats_txt(out, avg_vectors, avg_disp, se, se_mag, n_samples=10)
    lines = [l for l in out.read_text().splitlines() if l.strip()]
    # 4 header lines + n_ca data lines
    assert len(lines) == 4 + n_ca
