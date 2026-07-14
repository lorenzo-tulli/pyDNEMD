"""
GROMACS equilibration pipeline as a class.
"""
from pathlib import Path
from dnemd.gromacs import grompp, mdrun, genrestr, sed_posres
from dnemd.utils import ensure_dir, get_logger, copy_file

logger = get_logger("equilibration")

# stage -> (restraint group, force constant) or None for no restraint.
# "heavy" reuses the pdb2gmx-generated posre.itp (all protein heavy atoms);
# any other group is built via genrestr against the current index.ndx.
RESTRAINT_SCHEDULE = {
    "nvt1": ("heavy",    10),
    "nvt2": ("heavy",    5),
    "npt1": ("Backbone", 2),
    "npt2": None,
}
EQ_STAGES = list(RESTRAINT_SCHEDULE)


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
        """Run the full pipeline: EM -> nvt1-npt2 (restrained) -> production."""
        logger.info(f"=== Run {self.run_id} -> {self.root} ===")
        self.run_em()
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
        for itp in Path(self.cfg.topology).parent.glob("*.itp"):
            copy_file(itp, self.em_dir / itp.name)

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

    def _build_stage_restraint(self, stage_dir: Path, prev_gro: Path, group: str, fc: int):
        """
        Write posre.itp into stage_dir, restraining `group` at force constant `fc`.

        The topology's #ifdef POSRES block always includes a file literally
        named posre.itp, so each restrained stage needs its own copy under
        that exact name — reusing one shared file across stages is what
        silently pinned every stage to the same unscaled restraint before.
        """
        if group == "heavy":
            source = self.em_dir / "posre.itp"
            if not source.exists():
                logger.warning("posre.itp not found — heavy-atom restraints skipped.")
                return
            sed_posres(
                itp_in=str(source),
                itp_out=str(stage_dir / "posre.itp"),
                fc_value=str(fc),
            )
            return

        raw_itp = stage_dir / "posre_raw.itp"
        genrestr(
            gmx=self.cfg.gmx,
            gro=str(prev_gro),
            ndx="index.ndx",
            out_itp=raw_itp.name,
            group=group,
            cwd=stage_dir,
        )
        sed_posres(
            itp_in=str(raw_itp),
            itp_out=str(stage_dir / "posre.itp"),
            fc_value=str(fc),
        )

    def run_eq_stages(self):
        """Run nvt1 through npt2 sequentially, with per-stage position restraints."""
        for i, stage in enumerate(EQ_STAGES):
            stage_dir = ensure_dir(self.root / stage)
            self._copy_stage_inputs(stage_dir)

            prev_gro = (
                self.em_dir / "em.gro"
                if i == 0
                else self.root / EQ_STAGES[i - 1] / f"{EQ_STAGES[i - 1]}.gro"
            )

            restraint = RESTRAINT_SCHEDULE[stage]
            if restraint:
                group, fc = restraint
                self._build_stage_restraint(stage_dir, prev_gro, group, fc)

            grompp(
                gmx=self.cfg.gmx,
                mdp=str(self.mdp_dir / f"{stage}.mdp"),
                gro=str(prev_gro),
                top="topol.top",
                out_tpr=f"{stage}.tpr",
                ref_gro=str(prev_gro),
                ndx="index.ndx",
                cwd=stage_dir,
            )
            if not self.setup_only:
                mdrun(self.cfg.gmx, stage, cwd=stage_dir)

    def run_production(self):
        """Set up and run the unrestrained production simulation."""
        prod_dir   = ensure_dir(self.root / "prod")
        last_stage = EQ_STAGES[-1]
        last_gro   = self.root / last_stage / f"{last_stage}.gro"
        self._copy_stage_inputs(prod_dir)

        grompp(
            gmx=self.cfg.gmx,
            mdp=str(self.mdp_dir / "production.mdp"),
            gro=str(last_gro),
            top="topol.top",
            out_tpr="prod.tpr",
            ref_gro=str(last_gro),
            ndx="index.ndx",
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
