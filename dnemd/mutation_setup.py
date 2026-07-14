"""
Create NE/NP simulation input files for the mutation perturbation mode.
"""
from pathlib import Path
from dnemd.gromacs import grompp
from dnemd.utils import ensure_dir, get_logger, copy_file, run_piped

logger = get_logger("mutation_setup")


class MutationSetup:
    """
    Handles everything needed to create NE and NP input files for the
    mutation perturbation, from the hybrid topology (dnemd-create-hybrid-
    topology) and equilibrium trajectories.

    Same public interface as NESetup (check_required_files, create_inputs)
    so cli/create_NE_NP_files.py can dispatch on cfg.perturbation without
    needing separate code paths above this class.

    Unlike ligand removal, NE and NP both use the same hybrid topology
    throughout — only the MDP's lambda/velocity settings differ, so there's
    no separate "perturbed topology" concept to build here.

    The NE leg is a two-phase switch (see docs/design/mutation-perturbation.md):
    phase 1 ramps lambda 0->1, phase 2 continues at lambda=1 fixed. Only
    phase 1's input is built here — phase 2 needs phase 1's checkpoint,
    which doesn't exist until phase 1 has actually run, so building phase
    2's .tpr happens in run_ne_np.py at run time, not here at setup time.
    """

    def __init__(self, cfg):
        self.cfg = cfg

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def check_required_files(self) -> bool:
        missing = False
        for attr, label in [
            ("topology",  "topology"),
            ("input_gro", "input_gro"),
            ("index_ndx", "index_ndx"),
        ]:
            p = Path(getattr(self.cfg, attr))
            if not p.exists():
                logger.error(
                    f"Missing: {p} (cfg.{attr}) — run dnemd-create-hybrid-topology first?"
                )
                missing = True
        return not missing

    def create_inputs(self, run_ids: list[int], time_points_ns: list[int]):
        """Create NE (phase 1 only) and NP input files for all run/time combinations."""
        total = len(run_ids) * len(time_points_ns)
        done = 0
        for run_id in run_ids:
            for time_ns in time_points_ns:
                try:
                    self._create_ne_switch_input(run_id, time_ns)
                except Exception as e:
                    logger.warning(f"NE run {run_id}, {time_ns} ns failed: {e}")
                try:
                    self._create_np_input(run_id, time_ns)
                except Exception as e:
                    logger.warning(f"NP run {run_id}, {time_ns} ns failed: {e}")
                done += 1
                logger.info(f"Progress: {done}/{total}")

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _extract_frame(self, run_id: int, time_ns: int, out_gro: Path):
        """Dump a single frame from the EQ production trajectory."""
        prod_dir = Path(self.cfg.output_dir) / f"EQ_{run_id}" / "prod"
        run_piped(
            [self.cfg.gmx, "trjconv",
             "-f", str(prod_dir / "prod.xtc"),
             "-s", str(prod_dir / "prod.tpr"),
             "-o", str(out_gro),
             "-pbc", "whole",
             "-n", str(self.cfg.index_ndx),
             "-dump", str(time_ns),
             "-tu", "ns"],
            stdin_text="System\n",
            cwd=prod_dir,
        )

    def _copy_itps(self, leg_dir: Path):
        for itp in Path(self.cfg.topology).parent.glob("*.itp"):
            copy_file(itp, leg_dir / itp.name)

    def _create_ne_switch_input(self, run_id: int, time_ns: int):
        """
        Build phase 1 (switch) of the NE leg. Phase 2's .mdp is copied in
        now too (so it's on disk when run_ne_np.py needs it) but its .tpr
        isn't built until phase 1 has actually run.
        """
        leg_dir = ensure_dir(Path(self.cfg.output_dir) / f"NE_{run_id}" / f"{time_ns}ns")
        copy_file(Path(self.cfg.mdp_dir) / "switch_ph1.mdp", leg_dir / "switch_ph1.mdp")
        copy_file(Path(self.cfg.mdp_dir) / "switch_ph2.mdp", leg_dir / "switch_ph2.mdp")
        self._copy_itps(leg_dir)

        leg_gro = leg_dir / f"{time_ns}ns_NE.gro"
        logger.info(f"NE run {run_id}, {time_ns} ns: extracting frame...")
        self._extract_frame(run_id, time_ns, leg_gro)

        logger.info(f"NE run {run_id}, {time_ns} ns: building switch-phase input...")
        grompp(
            gmx=self.cfg.gmx,
            mdp="switch_ph1.mdp",
            gro=str(leg_gro),
            top=str(Path(self.cfg.topology)),
            out_tpr="MD_NE_switch.tpr",
            ref_gro=str(leg_gro),
            ndx=str(self.cfg.index_ndx),
            maxwarn=2,
            cwd=leg_dir,
        )

    def _create_np_input(self, run_id: int, time_ns: int):
        leg_dir = ensure_dir(Path(self.cfg.output_dir) / f"NP_{run_id}" / f"{time_ns}ns")
        copy_file(Path(self.cfg.mdp_dir) / "NP.mdp", leg_dir / "NP.mdp")
        self._copy_itps(leg_dir)

        leg_gro = leg_dir / f"{time_ns}ns_NP.gro"
        logger.info(f"NP run {run_id}, {time_ns} ns: extracting frame...")
        self._extract_frame(run_id, time_ns, leg_gro)

        logger.info(f"NP run {run_id}, {time_ns} ns: running grompp...")
        grompp(
            gmx=self.cfg.gmx,
            mdp="NP.mdp",
            gro=str(leg_gro),
            top=str(Path(self.cfg.topology)),
            out_tpr="MD_NP.tpr",
            ref_gro=str(leg_gro),
            ndx=str(self.cfg.index_ndx),
            maxwarn=2,
            cwd=leg_dir,
        )
