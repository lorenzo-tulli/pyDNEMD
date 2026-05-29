"""
Thin wrappers around GROMACS commands used by the pipeline.
All functions take an explicit gmx executable string so the caller
controls whether to use gmx or gmx_mpi.
"""
from pathlib import Path
from modules.utils import run, run_piped, ensure_dir, copy_file, get_logger

logger = get_logger("gromacs")


# ---------------------------------------------------------------------------
# Energy minimisation
# ---------------------------------------------------------------------------

def grompp(gmx: str, mdp: str, gro: str, top: str, out_tpr: str,
           ref_gro: str = None, maxwarn: int = 0, cwd: Path = None):
    cmd = [gmx, "grompp",
           "-f", mdp,
           "-c", gro,
           "-p", top,
           "-o", out_tpr,
           "-maxwarn", str(maxwarn)]
    if ref_gro:
        cmd += ["-r", ref_gro]
    run(cmd, cwd=cwd)


def mdrun(gmx: str, deffnm: str, cwd: Path = None, extra_flags: list = None):
    cmd = [gmx, "mdrun", "-v", "-deffnm", deffnm]
    if extra_flags:
        cmd += extra_flags
    run(cmd, cwd=cwd)


# ---------------------------------------------------------------------------
# Index / restraints
# ---------------------------------------------------------------------------

def make_index(gmx: str, gro: str, ndx_in: str, ndx_out: str,
               selection: str, cwd: Path = None):
    """Create a new index file with an added group defined by selection."""
    run_piped(
        [gmx, "make_ndx", "-f", gro, "-n", ndx_in, "-o", ndx_out],
        stdin_text=f"{selection}\nq\n",
        cwd=cwd,
    )


def genrestr(gmx: str, gro: str, ndx: str, out_itp: str,
             group: str = "CA", cwd: Path = None):
    """Generate position-restraint ITP for a specific atom group."""
    run_piped(
        [gmx, "genrestr", "-f", gro, "-n", ndx, "-o", out_itp],
        stdin_text=f"{group}\n",
        cwd=cwd,
    )


def sed_posres(itp_in: str, itp_out: str, fc_value: str):
    """Replace 1000 1000 1000 with fc_value in a posre ITP file."""
    text = Path(itp_in).read_text()
    text = text.replace("1000  1000  1000", f"{fc_value} {fc_value} {fc_value}")
    Path(itp_out).write_text(text)


# ---------------------------------------------------------------------------
# Trajectory tools
# ---------------------------------------------------------------------------

def trjconv_pbc(gmx: str, gro_in: str, tpr: str, gro_out: str,
                ndx: str, group: str, cwd: Path = None):
    """Correct PBC and write a new GRO file."""
    run_piped(
        [gmx, "trjconv",
         "-f", gro_in, "-s", tpr,
         "-o", gro_out,
         "-pbc", "whole",
         "-n", ndx],
        stdin_text=f"{group}\n",
        cwd=cwd,
    )


def trjconv_xtc(gmx: str, xtc_in: str, tpr: str, xtc_out: str,
                ndx: str, pbc_group: str, center_group: str,
                output_group: str, cwd: Path = None):
    """Center, correct PBC, and write XTC for analysis."""
    run_piped(
        [gmx, "trjconv",
         "-f", xtc_in, "-s", tpr,
         "-o", xtc_out,
         "-center", "-pbc", "mol", "-ur", "compact",
         "-n", ndx],
        stdin_text=f"{center_group}\n{output_group}\n",
        cwd=cwd,
    )


# ---------------------------------------------------------------------------
# RMSD / RMSF
# ---------------------------------------------------------------------------

def rms(gmx: str, xtc: str, tpr: str, out_xvg: str,
        ref_group: str, fit_group: str, ndx: str = None,
        cwd: Path = None):
    cmd = [gmx, "rms",
           "-f", xtc, "-s", tpr,
           "-o", out_xvg,
           "-tu", "ns"]
    if ndx:
        cmd += ["-n", ndx]
    run_piped(cmd, stdin_text=f"{ref_group}\n{fit_group}\n", cwd=cwd)


def rmsf(gmx: str, xtc: str, tpr: str, out_xvg: str,
         group: str, ndx: str = None, res: bool = True,
         cwd: Path = None):
    cmd = [gmx, "rmsf",
           "-f", xtc, "-s", tpr,
           "-o", out_xvg]
    if res:
        cmd += ["-res"]
    if ndx:
        cmd += ["-n", ndx]
    run_piped(cmd, stdin_text=f"{group}\n", cwd=cwd)
