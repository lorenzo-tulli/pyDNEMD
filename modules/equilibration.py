"""
GROMACS equilibration pipeline as a class.
"""
from pathlib import Path
from modules.gromacs import grompp, mdrun, genrestr, sed_posres
from modules.utils import ensure_dir, get_logger, copy_file

logger = get_logger("equilibration")

POSRES_FC = {"step1": 10, "step2": 5, "step3": 2, "step4": 1}
EQ_STAGES  = ["step1", "step2", "step3", "step4"]


class EquilibrationPipeline:
    """
    Manages one independent equilibration replicate.

    Parameters
    ----------
    cfg      : Config object
    run_id   : integer replicate number (e.g. 1, 2, 3)
    setup_only : if True, write input files but do not call mdrun
    """

    def __init__(self, cfg, run_id: int, setup_only: bool = False):
        self.cfg        = cfg
        self.run_id     = run_id
        self.setup_only = setup_only
        self.root       = ensure_dir(Path(cfg.output_dir) / f"EQ_{run_id}")
        self.mdp_dir    = Path(cfg.mdp_dir)
        self.em_dir     = ensure_dir(self.root / "em")

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def run_all(self):
        """Run the full pipeline: EM -> restraints -> step1-4 -> production."""
        logger.info(f"=== Run {self.run_id} -> {self.root} ===")
        self.run_em()
        self.build_restraints()
        self.run_eq_stages()
        self.run_production()
        logger.info(f"Run {self.run_id} {'setup' if self.setup_only else 'complete'}.")

    # ------------------------------------------------------------------
    # Pipeline stages
    # ------------------------------------------------------------------

    def run_em(self):
        """Copy inputs and run energy minimisation."""
        copy_file(self.cfg.input_gro, self.em_dir / "solv_ions.gro")
        copy_file(self.cfg.topology,  self.em_dir / "topol.top")
        copy_file(self.cfg.index_ndx, self.em_dir / "index.ndx")

        grompp(
            gmx=self.cfg.gmx,
            mdp=str(self.mdp_dir / "em.mdp"),
            gro="solv_ions.gro",
            top="topol.top",
            out_tpr="em.tpr",
            cwd=self.em_dir,
        )
        if not self.setup_only:
            mdrun(self.cfg.gmx, "em", cwd=self.em_dir)

    def build_restraints(self):
        """Generate posre_heavy.itp and posre_CA.itp."""
        # Heavy-atom restraints
        posre_src = self.em_dir / "posre.itp"
        if not posre_src.exists():
            posre_src = Path(self.cfg.topology).parent / "posre.itp"

        if posre_src.exists():
            sed_posres(
                itp_in=str(posre_src),
                itp_out=str(self.em_dir / "posre_heavy.itp"),
                fc_value=str(POSRES_FC["step1"]),
            )
        else:
            logger.warning("posre.itp not found — heavy-atom restraints skipped.")

        # Cα restraints
        gro = "em.gro" if not self.setup_only else "solv_ions.gro"
        genrestr(
            gmx=self.cfg.gmx,
            gro=gro,
            ndx="index.ndx",
            out_itp="posre_CA_raw.itp",
            group="CA",
            cwd=self.em_dir,
        )
        sed_posres(
            itp_in=str(self.em_dir / "posre_CA_raw.itp"),
            itp_out=str(self.em_dir / "posre_CA.itp"),
            fc_value=str(POSRES_FC["step2"]),
        )

    def run_eq_stages(self):
        """Run step1 through step4 sequentially."""
        for i, stage in enumerate(EQ_STAGES):
            stage_dir = ensure_dir(self.root / stage)
            self._copy_stage_inputs(stage_dir)

            prev_gro = (
                self.em_dir / "em.gro"
                if i == 0
                else self.root / EQ_STAGES[i - 1] / f"{EQ_STAGES[i - 1]}.gro"
            )

            grompp(
                gmx=self.cfg.gmx,
                mdp=str(self.mdp_dir / f"{stage}.mdp"),
                gro=str(prev_gro),
                top="topol.top",
                out_tpr=f"{stage}.tpr",
                ref_gro=str(prev_gro),
                cwd=stage_dir,
            )
            if not self.setup_only:
                mdrun(self.cfg.gmx, stage, cwd=stage_dir)

    def run_production(self):
        """Set up and run the unrestrained production simulation."""
        prod_dir  = ensure_dir(self.root / "prod")
        step4_gro = self.root / "step4" / "step4.gro"
        self._copy_stage_inputs(prod_dir)

        grompp(
            gmx=self.cfg.gmx,
            mdp=str(self.mdp_dir / "prod.mdp"),
            gro=str(step4_gro),
            top="topol.top",
            out_tpr="prod.tpr",
            ref_gro=str(step4_gro),
            cwd=prod_dir,
        )
        if not self.setup_only:
            mdrun(self.cfg.gmx, "prod", cwd=prod_dir)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _copy_stage_inputs(self, stage_dir: Path):
        """Copy all ITP files, topology, and index into a stage directory."""
        for itp in self.em_dir.glob("*.itp"):
            copy_file(itp, stage_dir / itp.name)
        copy_file(self.em_dir / "topol.top", stage_dir / "topol.top")
        copy_file(self.em_dir / "index.ndx", stage_dir / "index.ndx")
