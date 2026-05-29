import csv
import os
import tempfile
import numpy as np
import pytest
from dnemd.analysis import parse_xvg, write_summary_csv


XVG_CONTENT = """\
# This is a comment
@ title "RMSD"
@ xaxis label "Time (ns)"
@ yaxis label "RMSD (nm)"
0.0    0.10
1.0    0.20
2.0    0.15
3.0    0.18
"""


def test_parse_xvg_values(tmp_path):
    xvg = tmp_path / "test.xvg"
    xvg.write_text(XVG_CONTENT)
    x, y = parse_xvg(xvg)
    assert len(x) == 4
    assert len(y) == 4
    np.testing.assert_allclose(x, [0.0, 1.0, 2.0, 3.0])
    np.testing.assert_allclose(y, [0.10, 0.20, 0.15, 0.18])


def test_parse_xvg_skips_headers(tmp_path):
    xvg = tmp_path / "test.xvg"
    xvg.write_text("# comment\n@ header\n1.0  2.0\n")
    x, y = parse_xvg(xvg)
    assert len(x) == 1
    np.testing.assert_allclose(x[0], 1.0)
    np.testing.assert_allclose(y[0], 2.0)


def test_parse_xvg_empty(tmp_path):
    xvg = tmp_path / "empty.xvg"
    xvg.write_text("# only comments\n@ title\n")
    x, y = parse_xvg(xvg)
    assert len(x) == 0
    assert len(y) == 0


def test_write_summary_csv(tmp_path):
    rows = [
        {"run": 1, "mean_rmsd_A": 1.23, "max_rmsd_A": 2.34,
         "mean_rmsf_A": 0.56, "max_rmsf_A": 1.11},
        {"run": 2, "mean_rmsd_A": 1.45, "max_rmsd_A": 2.67,
         "mean_rmsf_A": 0.61, "max_rmsf_A": 1.22},
    ]
    out_csv = tmp_path / "summary.csv"
    write_summary_csv(rows, out_csv)
    assert out_csv.exists()
    with open(out_csv) as f:
        reader = list(csv.DictReader(f))
    assert len(reader) == 2
    assert reader[0]["run"] == "1"
    assert float(reader[1]["mean_rmsd_A"]) == pytest.approx(1.45)
