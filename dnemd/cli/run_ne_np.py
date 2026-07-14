#!/usr/bin/env python3
"""
Run a single NE or NP GROMACS simulation for one (run, time) combination.

The task ID maps to a (run_id, time_ns) pair using the same time points
defined by extract_start_ps / extract_frequency_ps / extract_end_ps in config.

With `perturbation: mutation`, the NE leg runs as two chained mdrun calls
(switch phase, then response phase built from the switch phase's
checkpoint) instead of one — see _run_mutation_ne(). The NP leg is
unaffected: single grompp+mdrun either way.

Usage (SLURM array):
    dnemd-run-ne --config config.yaml --task-id $SLURM_ARRAY_TASK_ID
    dnemd-run-np --config config.yaml --task-id $SLURM_ARRAY_TASK_ID

Usage (single run/time):
    dnemd-run-ne --config config.yaml --run 1 --time-ns 100
    dnemd-run-np --config config.yaml --run 1 --time-ns 100
"""
import argparse
import sys
from pathlib import Path

from dnemd.config import Config
from dnemd.gromacs import grompp, mdrun
from dnemd.utils import get_logger

logger = get_logger("run_ne_np")


def build_time_points(start_ps: int, frequency_ps: int, end_ps: int) -> list[int]:
    times, t = [], start_ps
    while t <= end_ps:
        times.append(t // 1000)
        t += frequency_ps
    return times


def run_simulation(leg: str, cfg: Config, run_id: int, time_ns: int):
    tpr_dir = Path(cfg.output_dir) / f"{leg}_{run_id}" / f"{time_ns}ns"

    if cfg.perturbation == "mutation" and leg == "NE":
        _run_mutation_ne(cfg, tpr_dir)
        return

    tpr = tpr_dir / f"MD_{leg}.tpr"
    if not tpr.exists():
        logger.error(f"TPR not found: {tpr}")
        sys.exit(1)

    logger.info(f"Running {leg} | run {run_id} | {time_ns} ns -> {tpr_dir}")
    mdrun(cfg.gmx, f"MD_{leg}", cwd=tpr_dir)


def _run_mutation_ne(cfg: Config, tpr_dir: Path):
    """
    Two-phase mutation NE: run the switch phase, then build the response
    phase's .tpr from the switch phase's checkpoint — this can't be built
    ahead of time in dnemd-create-ne-np, since the response phase needs
    the switch phase's final state, which doesn't exist until the switch
    phase has actually run. The response phase is deliberately named
    MD_NE (not MD_NE_ph2 or similar) so dnemd-extract's hardcoded
    MD_{leg}.xtc/.tpr lookup finds it unmodified.
    """
    switch_tpr = tpr_dir / "MD_NE_switch.tpr"
    if not switch_tpr.exists():
        logger.error(f"Switch-phase TPR not found: {switch_tpr}")
        sys.exit(1)

    logger.info(f"Running NE switch phase -> {tpr_dir}")
    mdrun(cfg.gmx, "MD_NE_switch", cwd=tpr_dir)

    logger.info("Building NE response phase from the switch-phase checkpoint...")
    grompp(
        gmx=cfg.gmx,
        mdp="switch_ph2.mdp",
        gro="MD_NE_switch.gro",
        top=str(Path(cfg.topology)),
        out_tpr="MD_NE.tpr",
        ref_gro="MD_NE_switch.gro",
        ndx=str(cfg.index_ndx),
        cpt="MD_NE_switch.cpt",
        maxwarn=2,
        cwd=tpr_dir,
    )

    logger.info("Running NE response phase...")
    mdrun(cfg.gmx, "MD_NE", cwd=tpr_dir)


def make_parser(leg: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=f"Run GROMACS {leg} simulation for one (run, time) combination."
    )
    parser.add_argument("--config",    required=True, help="Path to config.yaml")
    parser.add_argument(
        "--task-id", type=int, default=None,
        help="SLURM array task ID — mapped to (run, time_ns) automatically",
    )
    parser.add_argument("--run",     type=int, default=None,
                        help="Run ID (alternative to --task-id)")
    parser.add_argument("--time-ns", type=int, default=None,
                        help="Time point in ns (alternative to --task-id)")
    parser.add_argument("--start",     type=int, default=None,
                        help="First frame in ps (default: extract_start_ps from config)")
    parser.add_argument("--frequency", type=int, default=None,
                        help="Frame interval in ps (default: extract_frequency_ps from config)")
    parser.add_argument("--end",       type=int, default=None,
                        help="Last frame in ps (default: extract_end_ps from config)")
    return parser


def main(leg: str):
    parser = make_parser(leg)
    args = parser.parse_args()

    cfg = Config.from_yaml(args.config)

    start_ps     = args.start     or cfg.extract_start_ps
    frequency_ps = args.frequency or cfg.extract_frequency_ps
    end_ps       = args.end       or cfg.extract_end_ps

    if args.task_id is not None:
        time_points = build_time_points(start_ps, frequency_ps, end_ps)
        n_times     = len(time_points)
        run_id      = args.task_id // n_times + 1
        time_ns     = time_points[args.task_id % n_times]
    elif args.run is not None and args.time_ns is not None:
        run_id  = args.run
        time_ns = args.time_ns
    else:
        logger.error("Specify --task-id or both --run and --time-ns.")
        sys.exit(1)

    logger.info(f"Config       : {args.config}")
    logger.info(f"Leg          : {leg}")
    logger.info(f"Run ID       : {run_id}")
    logger.info(f"Time point   : {time_ns} ns")

    run_simulation(leg, cfg, run_id, time_ns)


def main_ne():
    main("NE")


def main_np():
    main("NP")