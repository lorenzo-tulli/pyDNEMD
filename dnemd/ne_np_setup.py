"""
Create NE/NP simulation input files.
"""
from pathlib import Path
from dnemd.gromacs import grompp
from dnemd.utils import ensure_dir, get_logger, copy_file, run_piped

logger = get_logger("ne_np_setup")


class NESetup:
    """
    Handles everything needed to create NE and NP input files from
    equilibrium trajectories.

    Parameters
    ----------
    cfg        : Config object
    perturb_dir: directory where the perturbed topology lives
    """

    def __init__(self, cfg, perturb_dir: Path = None):
        self.cfg            = cfg
        self.perturb_dir    = Path(perturb_dir or Path(cfg.output_dir) / "perturbed_topology")
        # User-provided: index.ndx with a Protein_Water_and_ions group added
        self.extraction_ndx = self.perturb_dir / "extraction_index.ndx"
        # Auto-generated: index built from perturbed_system.gro (no ligand)
        self.perturbed_ndx  = self.perturb_dir / "perturbed_index.ndx"

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def check_required_files(self) -> bool:
        """
        Check that extraction_index.ndx and topolperturb.top exist.
        Prints instructions and returns False if either is missing.
        """
        missing = False

        if not self.extraction_ndx.exists():
            logger.error(
                f"\nMissing: {self.extraction_ndx}\n\n"
                "Create an index file based on your production index.ndx, adding a group\n"
                "called 'Protein_Water_and_ions' that contains all atoms except the ligand.\n"
                "Then place it at the path shown above.\n\n"
                "Example (adjust group numbers to match your system):\n"
                f"  mkdir -p {self.perturb_dir}\n"
                f"  gmx_mpi make_ndx -f <output_dir>/EQ_1/em/em.gro \\\n"
                f"                   -n <input_dir>/index.ndx \\\n"
                f"                   -o {self.extraction_ndx}\n"
                "  # In the make_ndx prompt, list groups with Enter,\n"
                "  # then combine Protein and Water_and_ions by their group numbers, e.g.:\n"
                "  #   1 | 20\n"
                "  #   name <N> Protein_Water_and_ions\n"
                "  #   q\n"
            )
            missing = True

        topolperturb = self.perturb_dir / "topolperturb.top"
        if not topolperturb.exists():
            logger.error(
                f"\nMissing: {topolperturb}\n\n"
                "Create it by copying your topology file and removing the ligand\n"
                "from the [ molecules ] section, then place it at the path above.\n"
            )
            missing = True

        return not missing

    def build_perturbed_index(self):
        """
        Auto-generate perturbed_index.ndx from perturbed_system.gro.

        Steps:
          1. Extract Protein_Water_and_ions atoms from EQ_1/em/em.gro using
             extraction_index.ndx -> perturbed_system.gro
          2. Run make_ndx on perturbed_system.gro -> perturbed_index.ndx

        Skipped if perturbed_index.ndx already exists.
        """
        if self.perturbed_ndx.exists():
            logger.info("perturbed_index.ndx already exists — skipping rebuild.")
            return

        ensure_dir(self.perturb_dir)
        em_dir        = Path(self.cfg.output_dir) / "EQ_1" / "em"
        perturbed_gro = self.perturb_dir / "perturbed_system.gro"

        logger.info("Extracting perturbed system from EQ_1/em/em.gro...")
        run_piped(
            [self.cfg.gmx, "trjconv",
             "-f", str(em_dir / "em.gro"),
             "-s", str(em_dir / "em.tpr"),
             "-n", str(self.extraction_ndx),
             "-o", str(perturbed_gro)],
            stdin_text="Protein_Water_and_ions\n",
            cwd=self.perturb_dir,
        )

        logger.info("Building perturbed_index.ndx from perturbed_system.gro...")
        run_piped(
            [self.cfg.gmx, "make_ndx",
             "-f", str(perturbed_gro),
             "-o", str(self.perturbed_ndx)],
            stdin_text="q\n",
            cwd=self.perturb_dir,
        )
        logger.info(f"perturbed_index.ndx written to {self.perturbed_ndx}")

    def test_topology(self):
        """Validate the perturbed topology with grompp."""
        topolperturb = self.perturb_dir / "topolperturb.top"
        if not topolperturb.exists():
            logger.warning("topolperturb.top not found — skipping test.")
            return False

        grompp(
            gmx=self.cfg.gmx,
            mdp=str(Path(self.cfg.mdp_dir) / "Prod_RunNE.mdp"),
            gro=str(self.perturb_dir / "perturbed_system.gro"),
            top=str(topolperturb),
            out_tpr=str(self.perturb_dir / "test.tpr"),
            ref_gro=str(self.perturb_dir / "perturbed_system.gro"),
            ndx=str(self.perturbed_ndx),
            maxwarn=2,
            cwd=self.perturb_dir,
        )
        logger.info("Perturbed topology test passed.")
        return True

    def create_inputs(self, run_ids: list[int], time_points_ns: list[int]):
        """Create NE and NP input files for all run/time combinations."""
        self.build_perturbed_index()

        total = len(run_ids) * len(time_points_ns)
        done  = 0
        for run_id in run_ids:
            for time_ns in time_points_ns:
                for leg in ("NE", "NP"):
                    try:
                        self._create_leg_input(leg, run_id, time_ns)
                    except Exception as e:
                        logger.warning(f"{leg} run {run_id}, {time_ns} ns failed: {e}")
                done += 1
                logger.info(f"Progress: {done}/{total}")

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _extract_frame(self, run_id: int, time_ns: int, out_gro: Path,
                       ndx: str = None, group: str = "Protein_Water_and_ions"):
        """Dump a single frame from the production trajectory."""
        prod_dir = Path(self.cfg.output_dir) / f"EQ_{run_id}" / "prod"
        cmd = [self.cfg.gmx, "trjconv",
               "-f", str(prod_dir / "prod.xtc"),
               "-s", str(prod_dir / "prod.tpr"),
               "-o", str(out_gro),
               "-pbc", "whole",
               "-dump", str(time_ns),
               "-tu", "ns"]
        if ndx:
            cmd += ["-n", ndx]
        run_piped(cmd, stdin_text=f"{group}\n", cwd=prod_dir)

    def _create_leg_input(self, leg: str, run_id: int, time_ns: int):
        """
        Create input files for one leg (NE or NP) at one time point.

        NE: perturbed topology (ligand removed), keeps velocities from
            extracted frame (gen_vel = no in Prod_RunNE.mdp).
            - Frame extracted using extraction_index.ndx (group Protein_Water_and_ions)
            - grompp uses perturbed_index.ndx (built from perturbed structure)
        NP: original topology (ligand present), velocities reassigned from
            Maxwell-Boltzmann (gen_vel = yes in Prod_RunNP.mdp).
        """
        leg_dir = ensure_dir(
            Path(self.cfg.output_dir) / f"{leg}_{run_id}" / f"{time_ns}ns"
        )
        mdp_name = f"Prod_Run{leg}.mdp"
        copy_file(Path(self.cfg.mdp_dir) / mdp_name, leg_dir / mdp_name)

        leg_gro = leg_dir / f"{time_ns}ns_{leg}.gro"

        if leg == "NE":
            logger.info(f"NE run {run_id}, {time_ns} ns: extracting frame...")
            self._extract_frame(run_id, time_ns, leg_gro,
                                ndx=str(self.extraction_ndx),
                                group="Protein_Water_and_ions")
            top            = str(self.perturb_dir / "topolperturb.top")
            ndx_for_grompp = str(self.perturbed_ndx)
            for itp in self.perturb_dir.glob("*.itp"):
                copy_file(itp, leg_dir / itp.name)

        else:  # NP — unperturbed system, only velocities reassigned
            logger.info(f"NP run {run_id}, {time_ns} ns: extracting frame...")
            self._extract_frame(run_id, time_ns, leg_gro,
                                ndx=str(self.cfg.index_ndx),
                                group="System")
            top            = str(Path(self.cfg.topology))
            ndx_for_grompp = str(self.cfg.index_ndx)
            for itp in Path(self.cfg.topology).parent.glob("*.itp"):
                copy_file(itp, leg_dir / itp.name)

        tpr_name = f"MD_{leg}.tpr"
        logger.info(f"{leg} run {run_id}, {time_ns} ns: running grompp...")
        grompp(
            gmx=self.cfg.gmx,
            mdp=mdp_name,
            gro=str(leg_gro),
            top=top,
            out_tpr=tpr_name,
            ref_gro=str(leg_gro),
            ndx=ndx_for_grompp,
            maxwarn=2,
            cwd=leg_dir,
        )