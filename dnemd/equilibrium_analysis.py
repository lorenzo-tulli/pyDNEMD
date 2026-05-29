"""
Equilibrium simulation analysis: Cα RMSD and RMSF per run.
"""
from pathlib import Path
import numpy as np
from dnemd.gromacs import trjconv_xtc, rms, rmsf
from dnemd.analysis import parse_xvg
from dnemd.utils import ensure_dir, get_logger

logger = get_logger("equilibrium_analysis")


class EquilibriumAnalyser:
    """
    Computes Cα RMSD and RMSF for one equilibrium replicate.

    Parameters
    ----------
    cfg         : Config object
    run_id      : replicate number
    results_dir : where to write XVG and summary files
    """

    def __init__(self, cfg, run_id: int, results_dir: Path):
        self.cfg         = cfg
        self.run_id      = run_id
        self.results_dir = results_dir
        self.run_dir     = Path(cfg.output_dir) / f"EQ_{run_id}"
        self.prod_dir    = self.run_dir / "prod"
        self.em_dir      = self.run_dir / "em"

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def analyse(self) -> dict | None:
        """
        Run PBC correction, RMSD, and RMSF. Returns a summary dict,
        or None if required files are missing.
        """
        if not self._check_inputs():
            return None

        xtc_pbc = self._correct_pbc()
        self._compute_rmsd(xtc_pbc)
        self._compute_rmsf(xtc_pbc)
        return self._summarise()

    @property
    def rmsd_xvg(self) -> Path:
        return self.results_dir / f"rmsd_run{self.run_id}.xvg"

    @property
    def rmsf_xvg(self) -> Path:
        return self.results_dir / f"rmsf_run{self.run_id}.xvg"

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _check_inputs(self) -> bool:
        required = {
            "trajectory": self.prod_dir / "prod.xtc",
            "tpr":        self.prod_dir / "prod.tpr",
            "em.gro":     self.em_dir   / "em.gro",
        }
        for label, path in required.items():
            if not path.exists():
                logger.error(f"Missing {label} for run {self.run_id}: {path}")
                return False
        return True

    def _correct_pbc(self) -> Path:
        xtc_pbc = self.prod_dir / "prod_pbc.xtc"
        if not xtc_pbc.exists():
            logger.info(f"Run {self.run_id}: correcting PBC...")
            trjconv_xtc(
                gmx=self.cfg.gmx,
                xtc_in=str(self.prod_dir / "prod.xtc"),
                tpr=str(self.prod_dir / "prod.tpr"),
                xtc_out=str(xtc_pbc),
                ndx=str(self.em_dir / "index.ndx"),
                pbc_group="Protein",
                center_group="Protein",
                output_group="Protein",
                cwd=self.prod_dir,
            )
        else:
            logger.info(f"Run {self.run_id}: PBC-corrected trajectory already exists.")
        return xtc_pbc

    def _compute_rmsd(self, xtc_pbc: Path):
        logger.info(f"Run {self.run_id}: computing Cα RMSD...")
        rms(
            gmx=self.cfg.gmx,
            xtc=str(xtc_pbc),
            ref=str(self.em_dir / "em.gro"),
            out_xvg=str(self.rmsd_xvg),
            ref_group="C-alpha",
            fit_group="C-alpha",
            ndx=str(self.em_dir / "index.ndx"),
            cwd=self.prod_dir,
        )

    def _compute_rmsf(self, xtc_pbc: Path):
        logger.info(f"Run {self.run_id}: computing Cα RMSF...")
        rmsf(
            gmx=self.cfg.gmx,
            xtc=str(xtc_pbc),
            tpr=str(self.prod_dir / "prod.tpr"),
            out_xvg=str(self.rmsf_xvg),
            group="C-alpha",
            ndx=str(self.em_dir / "index.ndx"),
            res=True,
            cwd=self.prod_dir,
        )

    def _summarise(self) -> dict:
        _, r_rmsd = parse_xvg(self.rmsd_xvg)
        _, r_rmsf = parse_xvg(self.rmsf_xvg)
        r_rmsd_A  = r_rmsd * 10
        r_rmsf_A  = r_rmsf * 10
        return {
            "run":         self.run_id,
            "mean_rmsd_A": round(float(np.mean(r_rmsd_A)), 3),
            "max_rmsd_A":  round(float(np.max(r_rmsd_A)),  3),
            "mean_rmsf_A": round(float(np.mean(r_rmsf_A)), 3),
            "max_rmsf_A":  round(float(np.max(r_rmsf_A)),  3),
        }
