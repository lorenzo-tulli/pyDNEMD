"""
Create NE/NP simulation input files.
"""
import shutil
from pathlib import Path
from dnemd.gromacs import make_index, trjconv_pbc, grompp
from dnemd.utils import ensure_dir, get_logger, copy_file, run, run_piped

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
        self.cfg         = cfg
        self.perturb_dir = perturb_dir or Path(cfg.output_dir) / "perturbed_topology"
        self.ndx         = self.perturb_dir / "perturbed_index.ndx"

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def check_required_files(self) -> bool:
        """
        Check that perturbed_index.ndx and topolperturb.top exist.
        Prints instructions and returns False if either is missing.
        """
        missing = False

        if not self.ndx.exists():
            logger.error(
                f"\nMissing: {self.ndx}\n"
                "Create a custom index file with a group called 'Protein_Water_and_ions'\n"
                "containing all atoms except the ligand, then place it at the path above.\n\n"
                "Example (adjust group numbers for your system):\n"
                f"  mkdir -p {self.perturb_dir}\n"
                f"  gmx_mpi make_ndx -f <output_dir>/EQ_1/em/em.gro \\\n"
                f"                   -n <input_dir>/index.ndx \\\n"
                f"                   -o {self.ndx}\n"
                "  # In the make_ndx prompt, list groups with Enter, then combine\n"
                "  # Protein and Water_and_ions using their group numbers, e.g.:\n"
                "  #   1 | 20\n"
                "  #   q\n"
            )
            missing = True

        topolperturb = self.perturb_dir / "topolperturb.top"
        if not topolperturb.exists():
            logger.error(
                f"\nMissing: {topolperturb}\n"
                "Create it by copying your topology file and removing the ligand\n"
                "from the [ molecules ] section, then place it at the path above.\n"
            )
            missing = True

        return not missing

    def test_topology(self):
        """Validate the perturbed topology with grompp."""
        topolperturb = self.perturb_dir / "topolperturb.top"
        if not topolperturb.exists():
            logger.warning("topolperturb.top not found — skipping test.")
            return False

        grompp(
            gmx=self.cfg.gmx,
            mdp=str(Path(self.cfg.mdp_dir) / "Prod_RunNE.mdp"),
            gro=str(self.perturb_dir / "perturbed_em.gro"),
            top=str(topolperturb),
            out_tpr=str(self.perturb_dir / "test.tpr"),
            ref_gro=str(self.perturb_dir / "perturbed_em.gro"),
            maxwarn=2,
            cwd=self.perturb_dir,
        )
        logger.info("Perturbed topology test passed.")
        return True

    def create_inputs(self, run_ids: list[int], time_points_ns: list[int]):
        """Create NE and NP input files for all run/time combinations."""
        total = len(run_ids) * len(time_points_ns)
        done  = 0
        for run_id in run_ids:
            for time_ns in time_points_ns:
                try:
                    self._create_leg_input("NE", run_id, time_ns)
                    self._create_leg_input("NP", run_id, time_ns)
                except Exception as e:
                    logger.warning(f"Run {run_id}, {time_ns} ns failed: {e}")
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
                                ndx=str(self.ndx),
                                group="Protein_Water_and_ions")
            top = str(self.perturb_dir / "topolperturb.top")
            ndx_for_grompp = str(self.ndx)
            # Copy ITP files needed by topolperturb.top
            for itp in self.perturb_dir.glob("*.itp"):
                copy_file(itp, leg_dir / itp.name)

        else:  # NP — unperturbed system, only velocities reassigned
            logger.info(f"NP run {run_id}, {time_ns} ns: extracting frame...")
            self._extract_frame(run_id, time_ns, leg_gro,
                                ndx=str(self.cfg.index_ndx),
                                group="System")
            top = str(Path(self.cfg.topology))
            ndx_for_grompp = str(self.cfg.index_ndx)
            # Copy ITP files needed by topol.top
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