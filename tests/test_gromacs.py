import pytest
from dnemd.gromacs import sed_posres


def test_sed_posres_double_space(tmp_path):
    itp_in  = tmp_path / "in.itp"
    itp_out = tmp_path / "out.itp"
    itp_in.write_text("[ position_restraints ]\n1  1  1000  1000  1000\n")
    sed_posres(str(itp_in), str(itp_out), "500")
    content = itp_out.read_text()
    assert "1000  1000  1000" not in content
    assert "500" in content


def test_sed_posres_multiple_spaces(tmp_path):
    itp_in  = tmp_path / "in.itp"
    itp_out = tmp_path / "out.itp"
    itp_in.write_text("[ position_restraints ]\n1  1  1000    1000    1000\n")
    sed_posres(str(itp_in), str(itp_out), "200")
    content = itp_out.read_text()
    assert "1000    1000    1000" not in content
    assert "200" in content


def test_sed_posres_four_spaces(tmp_path):
    itp_in  = tmp_path / "in.itp"
    itp_out = tmp_path / "out.itp"
    itp_in.write_text("1  1  1000    1000    1000\n2  1  1000    1000    1000\n")
    sed_posres(str(itp_in), str(itp_out), "10")
    content = itp_out.read_text()
    assert content.count("1000") == 0


def test_sed_posres_no_match_leaves_file_unchanged(tmp_path):
    itp_in  = tmp_path / "in.itp"
    itp_out = tmp_path / "out.itp"
    original = "[ position_restraints ]\n1  1  500  500  500\n"
    itp_in.write_text(original)
    sed_posres(str(itp_in), str(itp_out), "200")
    assert itp_out.read_text() == original


def test_sed_posres_preserves_other_content(tmp_path):
    itp_in  = tmp_path / "in.itp"
    itp_out = tmp_path / "out.itp"
    itp_in.write_text(
        "[ position_restraints ]\n"
        "; atom  type  fx   fy   fz\n"
        "1  1  1000  1000  1000\n"
    )
    sed_posres(str(itp_in), str(itp_out), "5")
    content = itp_out.read_text()
    assert "[ position_restraints ]" in content
    assert "; atom  type  fx   fy   fz" in content
