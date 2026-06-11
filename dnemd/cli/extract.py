#!/usr/bin/env python3
"""
Extract EQ, NE, and NP frames needed for DNEMD analysis.

Must be run AFTER all equilibrium, NE, and NP simulations have finished
and BEFORE running analyse_dnemd.py.

This script mirrors trjdump.sh:
  - EQ: extracts a PDB at each ps time point inside the <ns>ns window,
         rot+trans fitted to the reference group
  - NE: centers -> fixes PBC -> dumps each ps frame with rot+trans fit
  - NP: same pipeline as NE

The ns windows to extract are submitted as a SLURM array (one task per ns
window) or run sequentially for all windows with --all.

Usage (SLURM array, one task per ns window):
    sbatch --array=0-39 scripts/extract.py --config examples/config.yaml

Usage (sequential, all ns windows):
    python scripts/extract.py --config examples/config.yaml --all

Usage (single ns window):
    python scripts/extract.py --config examples/config.yaml --ns 100

Required config.yaml additions (see examples/config.yaml):
    eq_sim_dir:  /path/to/equil_simuls
    ne_sim_dir:  /path/to/D-NEMD_simuls
    np_sim_dir:  /path/to/null-perturbations_simuls
    ne_index:    /path/to/ne_np_index.ndx
    extract_ns_start:    50      # ns — first window
    extract_ns_end:      250     # ns — last window
    extract_ns_interval: 5       # ns — step between windows
"""
import argparse
import os
import sys
from pathlib import Path

from dnemd.config import Config
from dnemd.trjdump import (
    TrajectoryDumper,
    DEFAULT_PS_TIMEPOINTS,
)
from dnemd.utils import get_logger

logger = get_logger("extract")


def get_ns_windows(cfg: Config) -> list[int]:
    """Build the list of ns windows from config."""
    start    = getattr(cfg, "extract_ns_start",    50)
    end      = getattr(cfg, "extract_ns_end",      250)
    interval = getattr(cfg, "extract_ns_interval", 5)
    return list(range(start, end + 1, interval))


def resolve_ns(args, cfg: Config) -> list[int]:
    """Determine which ns windows to process based on CLI flags."""
    all_windows = get_ns_windows(cfg)

    if args.ns is not None:
        return [args.ns]

    if args.all:
        return all_windows

    # SLURM array mode
    task_id = args.task_id
    if task_id is None:
        task_id_env = os.environ.get("SLURM_ARRAY_TASK_ID")
        if task_id_env is not None:
            task_id = int(task_id_env)

    if task_id is not None:
        if task_id >= len(all_windows):
            logger.error(
                f"Task ID {task_id} out of range "
                f"(there are {len(all_windows)} ns windows: "
                f"{all_windows[0]}–{all_windows[-1]} ns)."
            )
            sys.exit(1)
        return [all_windows[task_id]]

    # No mode selected
    logger.error(
        "Specify one of: --ns <value>, --all, --task-id <id>, "
        "or set SLURM_ARRAY_TASK_ID."
    )
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract EQ/NE/NP frames for DNEMD analysis. "
            "Run after all simulations finish, before analyse_dnemd.py."
        )
    )
    parser.add_argument("--config",   required=True, help="Path to config.yaml")
    parser.add_argument("--ns",       type=int, default=None,
                        help="Single ns window to extract (e.g. 100)")
    parser.add_argument("--all",      action="store_true",
                        help="Extract all ns windows sequentially")
    parser.add_argument("--task-id",  type=int, default=None,
                        help="SLURM array task ID (index into ns window list)")
    parser.add_argument(
        "--leg", choices=["eq", "ne", "np", "all"], default="all",
        help="Which leg to extract (default: all)"
    )
    parser.add_argument(
        "--ps-timepoints", nargs="+", type=int, default=None,
        help="ps time points to extract within each window (default: from config time_points_ps)",
    )
    parser.add_argument(
        "--runs", nargs="+", type=int, default=None,
        help="Run IDs to process (default: 1..n_runs from config)"
    )
    parser.add_argument(
        "--fit-group",    default="1",
        help="GROMACS index group number for rot+trans fitting (default: 1 = Protein)"
    )
    parser.add_argument(
        "--output-group", default="1",
        help="GROMACS index group number for output (default: 1 = Protein)"
    )
    parser.add_argument(
        "--center-group", default="1",
        help="GROMACS index group number for centering NE/NP (default: 1 = Protein)"
    )
    args = parser.parse_args()

    # ---------------------------------------------------------------- Config --
    cfg = Config.from_yaml(args.config)

    # Resolve required paths from config
    eq_sim_dir = Path(getattr(cfg, "eq_sim_dir", ""))
    ne_sim_dir = Path(getattr(cfg, "ne_sim_dir", ""))
    np_sim_dir = Path(getattr(cfg, "np_sim_dir", ""))
    ne_index   = Path(getattr(cfg, "ne_index",   ""))

    if args.leg in ("ne", "np", "all"):
        if not ne_index.exists():
            logger.error(
                f"NE/NP index file not found: {ne_index}\n"
                "Set 'ne_index' in config.yaml to the path of your index.ndx."
            )
            sys.exit(1)

    # Resolve run IDs and time points
    runs          = args.runs if args.runs else list(range(1, cfg.n_runs + 1))
    ps_timepoints = args.ps_timepoints or cfg.time_points_ps or DEFAULT_PS_TIMEPOINTS
    ns_windows   = resolve_ns(args, cfg)

    logger.info(f"ns windows    : {ns_windows}")
    logger.info(f"ps time points: {ps_timepoints}")
    logger.info(f"runs          : {runs}")
    logger.info(f"leg           : {args.leg}")

    # --------------------------------------------------------------- Extract --
    for ns in ns_windows:
        logger.info(f"======= ns window: {ns} ns =======")

        if args.leg in ("eq", "all") and eq_sim_dir.exists():
            TrajectoryDumper(cfg.gmx, "EQ", eq_sim_dir,
                            fit_group=args.fit_group,
                            output_group=args.output_group).dump(runs, ns, ps_timepoints)

        if args.leg in ("ne", "all") and ne_sim_dir.exists():
            TrajectoryDumper(cfg.gmx, "NE", ne_sim_dir, ne_index,
                            fit_group=args.fit_group,
                            output_group=args.output_group,
                            center_group=args.center_group).dump(runs, ns, ps_timepoints)

        if args.leg in ("np", "all") and np_sim_dir.exists():
            TrajectoryDumper(cfg.gmx, "NP", np_sim_dir, ne_index,
                            fit_group=args.fit_group,
                            output_group=args.output_group,
                            center_group=args.center_group).dump(runs, ns, ps_timepoints)


    logger.info("extract.py finished.")


if __name__ == "__main__":
    main()
