import pytest
from dnemd.config import Config


def test_default_values(tmp_path):
    cfg_file = tmp_path / "config.yaml"
    cfg_file.write_text("system_name: TEST\n")
    cfg = Config.from_yaml(str(cfg_file))
    assert cfg.system_name == "TEST"
    assert cfg.gmx == "gmx_mpi"
    assert cfg.n_runs == 5
    assert cfg.se_threshold == 1
    assert cfg.extract_ns_start == 50
    assert cfg.extract_ns_end == 250
    assert cfg.extract_ns_interval == 5


def test_override_values(tmp_path):
    cfg_file = tmp_path / "config.yaml"
    cfg_file.write_text(
        "system_name: MYPROTEIN\n"
        "gmx: gmx\n"
        "n_runs: 3\n"
        "se_threshold: 2\n"
    )
    cfg = Config.from_yaml(str(cfg_file))
    assert cfg.system_name == "MYPROTEIN"
    assert cfg.gmx == "gmx"
    assert cfg.n_runs == 3
    assert cfg.se_threshold == 2


def test_extraction_fields_loaded(tmp_path):
    cfg_file = tmp_path / "config.yaml"
    cfg_file.write_text(
        "eq_sim_dir: /path/to/eq\n"
        "ne_sim_dir: /path/to/ne\n"
        "np_sim_dir: /path/to/np\n"
        "ne_index: /path/to/index.ndx\n"
        "extract_ns_start: 100\n"
        "extract_ns_end: 300\n"
        "extract_ns_interval: 10\n"
    )
    cfg = Config.from_yaml(str(cfg_file))
    assert cfg.eq_sim_dir == "/path/to/eq"
    assert cfg.ne_sim_dir == "/path/to/ne"
    assert cfg.np_sim_dir == "/path/to/np"
    assert cfg.ne_index == "/path/to/index.ndx"
    assert cfg.extract_ns_start == 100
    assert cfg.extract_ns_end == 300
    assert cfg.extract_ns_interval == 10


def test_unknown_fields_ignored(tmp_path):
    cfg_file = tmp_path / "config.yaml"
    cfg_file.write_text("system_name: TEST\nunknown_field: 999\n")
    cfg = Config.from_yaml(str(cfg_file))
    assert cfg.system_name == "TEST"
    assert not hasattr(cfg, "unknown_field")


def test_time_points_ps(tmp_path):
    cfg_file = tmp_path / "config.yaml"
    cfg_file.write_text("time_points_ps: [0, 100, 1000]\n")
    cfg = Config.from_yaml(str(cfg_file))
    assert cfg.time_points_ps == [0, 100, 1000]
