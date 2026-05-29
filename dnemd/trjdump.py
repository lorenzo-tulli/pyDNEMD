"""
Frame extraction for DNEMD analysis.
"""
from pathlib import Path
from dnemd.utils import run_piped, ensure_dir, get_logger

logger = get_logger("trjdump")

DEFAULT_PS_TIMEPOINTS = [
    0, 100, 200, 300, 400, 500, 600, 700, 800, 900,
    1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000,
]


class TrajectoryDumper:
    """
    Extracts PDB frames from EQ, NE, or NP trajectories for DNEMD analysis.

    Parameters
    ----------
    gmx          : GROMACS executable name
    leg          : one of "EQ", "NE", "NP"
    sim_path     : root directory of the simulation leg
    index_ndx    : path to index file (required for NE/NP, ignored for EQ)
    fit_group    : GROMACS group number for rot+trans fitting
    output_group : GROMACS group number for output atoms
    center_group : GROMACS group number for centering (NE/NP only)
    """

    def __init__(
        self,
        gmx: str,
        leg: str,
        sim_path: Path,
        index_ndx: Path = None,
        fit_group: str = "1",
        output_group: str = "1",
        center_group: str = "1",
    ):
        if leg not in ("EQ", "NE", "NP"):
            raise ValueError(f"leg must be 'EQ', 'NE', or 'NP', got '{leg}'")
        self.gmx          = gmx
        self.leg          = leg
        self.sim_path     = Path(sim_path)
        self.index_ndx    = Path(index_ndx) if index_ndx else None
        self.fit_group    = fit_group
        self.output_group = output_group
        self.center_group = center_group
        self.trjdump_root = ensure_dir(sim_path / f"TRJDUMP_{leg}")

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def dump(self, runs: list[int], ns: int, ps_timepoints: list[int] = None):
        """Extract frames for all runs at the given ns window."""
        ps_timepoints = ps_timepoints or DEFAULT_PS_TIMEPOINTS
        for run in runs:
            self._process_run(run, ns, ps_timepoints)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _run_dir(self, run: int, ns: int) -> Path:
        return self.sim_path / f"r{run}" / f"{ns}ns"

    def _out_dir(self, run: int, ns: int) -> Path:
        return ensure_dir(
            self.trjdump_root / f"{self.leg}_{run}" / f"{ns}ns_{self.leg}"
        )

    def _locate_trajectories(self, run: int, ns: int):
        """Return (xtc, tpr) paths or raise a warning and return (None, None)."""
        run_dir = self._run_dir(run, ns)

        if self.leg == "EQ":
            xtc = run_dir / f"r{run}.xtc"
            tpr = run_dir / f"r{run}.tpr"
        else:
            xtc = run_dir / f"r{run}-{ns}ns.xtc"
            tpr = run_dir / "md.tpr"

        if not xtc.exists() or not tpr.exists():
            logger.warning(
                f"{self.leg} run {run}, {ns} ns: trajectory not found, skipping."
            )
            return None, None
        return xtc, tpr

    def _process_run(self, run: int, ns: int, ps_timepoints: list[int]):
        xtc, tpr = self._locate_trajectories(run, ns)
        if xtc is None:
            return

        out_dir = self._out_dir(run, ns)
        ndx     = str(self.index_ndx) if self.index_ndx else None

        if self.leg == "EQ":
            self._dump_eq_frames(run, ns, xtc, tpr, out_dir, ps_timepoints)
        else:
            xtc_pbc = self._prepare_trajectory(run, ns, xtc, tpr, ndx)
            if xtc_pbc:
                self._dump_ne_np_frames(run, ns, xtc_pbc, tpr, out_dir, ps_timepoints, ndx)

    def _prepare_trajectory(self, run: int, ns: int, xtc_raw: Path,
                             tpr: Path, ndx: str) -> Path | None:
        """Center and fix PBC. Returns path to the cleaned trajectory."""
        run_dir    = self._run_dir(run, ns)
        xtc_center = run_dir / "center.xtc"
        xtc_pbc    = run_dir / f"{self.leg}.xtc"

        if not xtc_center.exists():
            logger.info(f"{self.leg} run {run} | {ns} ns: centering...")
            run_piped(
                [self.gmx, "trjconv",
                 "-f", str(xtc_raw), "-s", str(tpr),
                 "-o", str(xtc_center), "-center", "-n", ndx],
                stdin_text=f"{self.center_group}\n{self.output_group}\n",
            )

        if not xtc_pbc.exists():
            logger.info(f"{self.leg} run {run} | {ns} ns: fixing PBC...")
            run_piped(
                [self.gmx, "trjconv",
                 "-f", str(xtc_center), "-s", str(tpr),
                 "-o", str(xtc_pbc), "-pbc", "mol", "-n", ndx],
                stdin_text=f"{self.output_group}\n",
            )

        return xtc_pbc if xtc_pbc.exists() else None

    def _dump_eq_frames(self, run: int, ns: int, xtc: Path, tpr: Path,
                         out_dir: Path, ps_timepoints: list[int]):
        for ps in ps_timepoints:
            out_pdb = out_dir / f"Run{run}EQ{ns}ns{ps}ps.pdb"
            if out_pdb.exists():
                logger.info(f"Already exists, skipping: {out_pdb.name}")
                continue

            start_time = (ns * 1000) - 1
            end_time   = start_time + ps + 1
            logger.info(f"EQ run {run} | {ns} ns | {ps} ps -> {out_pdb.name}")
            run_piped(
                [self.gmx, "trjconv",
                 "-f", str(xtc), "-s", str(tpr),
                 "-o", str(out_pdb),
                 "-b", str(start_time), "-dump", str(end_time),
                 "-fit", "rot+trans"],
                stdin_text=f"{self.fit_group}\n{self.output_group}\n",
            )

    def _dump_ne_np_frames(self, run: int, ns: int, xtc_pbc: Path, tpr: Path,
                            out_dir: Path, ps_timepoints: list[int], ndx: str):
        for ps in ps_timepoints:
            out_pdb = out_dir / f"Run{run}{self.leg}{ns}ns{ps}ps.pdb"
            if out_pdb.exists():
                logger.info(f"Already exists, skipping: {out_pdb.name}")
                continue

            logger.info(f"{self.leg} run {run} | {ns} ns | {ps} ps -> {out_pdb.name}")
            run_piped(
                [self.gmx, "trjconv",
                 "-f", str(xtc_pbc), "-s", str(tpr),
                 "-o", str(out_pdb),
                 "-dump", str(ps), "-fit", "rot+trans", "-n", ndx],
                stdin_text=f"{self.fit_group}\n{self.output_group}\n",
            )
