"""
Create NE/NP simulation input files.
"""
import shutil
from pathlib import Path
from modules.gromacs import make_index, trjconv_pbc, grompp
from modules.utils import ensure_dir, get_logger, copy_file, run, run_piped

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

    def build_topology(self, run_id: int = 1):
        """Create perturbed_em.gro and perturbed_index.ndx."""
        ensure_dir(self.perturb_dir)
        em_dir = Path(self.cfg.output_dir) / f"EQ_{run_id}" / "em"
        em_gro = em_dir / "em.gro"
        em_tpr = em_dir / "em.tpr"

        if not em_gro.exists():
            raise FileNotFoundError(
                f"EM structure not found: {em_gro}\n"
                "Run run_equilibrium.py first."
            )

        shutil.copy2(em_dir / "index.ndx", self.ndx)
        logger.info("Building Protein_Solvent index group...")
        make_index(
            gmx=self.cfg.gmx,
            gro=str(em_gro),
            ndx_in=str(self.ndx),
            ndx_out=str(self.ndx),
            selection="Protein | Water_and_ions",
            cwd=self.perturb_dir,
        )

        perturbed_gro = self.perturb_dir / "perturbed_em.gro"
        logger.info("Extracting protein+solvent structure...")
        trjconv_pbc(
            gmx=self.cfg.gmx,
            gro_in=str(em_gro),
            tpr=str(em_tpr),
            gro_out=str(perturbed_gro),
            ndx=str(self.ndx),
            group="Protein_Water_and_ions",
            cwd=self.perturb_dir,
        )
        logger.info(
            "ACTION REQUIRED: create topolperturb.top in "
            f"{self.perturb_dir} by removing ligand entries from topol.top."
        )

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

    def _extract_frame(self, run_id: int, time_ns: int, out_gro: Path):
        """Dump a single frame from the production trajectory."""
        prod_dir = Path(self.cfg.output_dir) / f"EQ_{run_id}" / "prod"
        run_piped(
            [self.cfg.gmx, "trjconv",
             "-f", str(prod_dir / "prod.xtc"),
             "-s", str(prod_dir / "prod.tpr"),
             "-o", str(out_gro),
             "-pbc", "whole",
             "-dump", str(time_ns),
             "-tu", "ns",
             "-n", str(self.ndx)],
            stdin_text="Protein_Water_and_ions\n",
            cwd=prod_dir,
        )

    def _create_leg_input(self, leg: str, run_id: int, time_ns: int):
        """
        Create input files for one leg (NE or NP) at one time point.
        NE and NP share the extracted GRO; NP uses a different MDP
        (gen_vel = yes) to reassign velocities.
        """
        leg_dir = ensure_dir(
            Path(self.cfg.output_dir) / f"TRJDUMP_{leg}"
            / f"{leg}_{run_id}" / f"{time_ns}ns"
        )
        mdp_name = f"Prod_Run{leg}.mdp"
        copy_file(Path(self.cfg.mdp_dir) / mdp_name, leg_dir / mdp_name)

        posre = self.perturb_dir / "posre.itp"
        if posre.exists():
            copy_file(posre, leg_dir / "posre.itp")

        # NE extracts its own frame; NP reuses it
        ne_gro = (
            Path(self.cfg.output_dir) / "TRJDUMP_NE"
            / f"NE_{run_id}" / f"{time_ns}ns" / f"{time_ns}ns_NE.gro"
        )
        leg_gro = leg_dir / f"{time_ns}ns_{leg}.gro"

        if leg == "NE":
            logger.info(f"NE run {run_id}, {time_ns} ns: extracting frame...")
            self._extract_frame(run_id, time_ns, leg_gro)
        else:
            if not ne_gro.exists():
                logger.warning(f"NE gro not found for NP: {ne_gro}. Skipping.")
                return
            shutil.copy2(ne_gro, leg_gro)

        tpr_name = f"MD_{leg}.tpr"
        logger.info(f"{leg} run {run_id}, {time_ns} ns: running grompp...")
        grompp(
            gmx=self.cfg.gmx,
            mdp=mdp_name,
            gro=str(leg_gro),
            top=str(self.perturb_dir / "topolperturb.top"),
            out_tpr=tpr_name,
            ref_gro=str(leg_gro),
            maxwarn=2,
            cwd=leg_dir,
        )