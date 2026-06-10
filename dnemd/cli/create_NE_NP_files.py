#!/usr/bin/env python3

"""
1. Auto-generate perturbed_index.ndx from the first EQ energy-minimised structure.
2. Extract frames from each equilibrium run at the requested interval.
3. Create GROMACS input files for Non-Equilibrium (NE) simulations
   (perturbed topology, velocities NOT reassigned).
4. Create GROMACS input files for Null-Perturbation (NP) simulations
   (original topology, velocities reassigned from Maxwell-Boltzmann distribution).

Files the user must place in output/perturbed_topology/ before running:

    extraction_index.ndx  — copy of your production index.ndx with an extra group
                             called 'Protein_Water_and_ions' (protein + water + ions,
                             ligand excluded).  Create it with:
                               gmx_mpi make_ndx -f output/EQ_1/em/em.gro \\
                                                -n inputs/index.ndx \\
                                                -o output/perturbed_topology/extraction_index.ndx
                             then in the make_ndx prompt combine the groups, e.g.:
                               1 | 20
                               name <N> Protein_Water_and_ions
                               q

    topolperturb.top      — copy of topol.top with the ligand removed from [ molecules ]:
                               Protein  1
                               SOL      39553
                               NA       119
                               CL       114

Directory layout created:
    output/
    ├── perturbed_topology/
    │   ├── extraction_index.ndx   ← user must provide
    │   ├── topolperturb.top       ← user must provide
    │   ├── perturbed_system.gro   ← auto-generated
    │   └── perturbed_index.ndx    ← auto-generated
    ├── NE_<run>/<time>ns/
    │   ├── <time>ns_NE.gro
    │   ├── Prod_RunNE.mdp
    │   └── MD_NE.tpr
    └── NP_<run>/<time>ns/
        ├── <time>ns_NP.gro
        ├── Prod_RunNP.mdp
        └── MD_NP.tpr

Usage:
    dnemd-create-ne-np --config config.yaml

    --start      : first frame to extract, in ps  (default: extract_start_ps from config)
    --frequency  : interval between frames, in ps (default: extract_frequency_ps from config)
    --end        : last frame in ps               (default: extract_end_ps from config)
"""

import argparse
from pathlib import Path
from dnemd.config import Config
from dnemd.ne_np_setup import NESetup
from dnemd.utils import get_logger, run

logger = get_logger("create_NE_NP_files")


def get_production_length_ps(tpr: str, gmx: str) -> int:
    try:
        result = run([gmx, "check", "-f", tpr], check=False)
        for line in result.stdout.splitlines():
            if "Last frame" in line or "nsteps" in line.lower():
                parts = line.split()
                for p in parts:
                    if p.isdigit():
                        return int(p)
    except Exception:
        pass
    logger.warning("Could not determine production length — defaulting to 500 ns.")
    return 500_000


def frame_times_ns(start_ps: int, frequency_ps: int, end_ps: int) -> list[int]:
    times, t = [], start_ps
    while t <= end_ps:
        times.append(t // 1000)
        t += frequency_ps
    return times


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config",                required=True)
    parser.add_argument("--start",     type=int,   default=None,
                        help="First frame to extract in ps (default: extract_start_ps from config)")
    parser.add_argument("--frequency", type=int,   default=None,
                        help="Interval between frames in ps (default: extract_frequency_ps from config)")
    parser.add_argument("--end",       type=int,   default=None,
                        help="Last frame in ps (default: extract_end_ps from config, or auto-detected from TPR)")
    parser.add_argument("--run",       type=int,   default=None)
    parser.add_argument("--skip-topology-test",  action="store_true")
    args = parser.parse_args()

    cfg     = Config.from_yaml(args.config)
    setup   = NESetup(cfg)
    run_ids = [args.run] if args.run else list(range(1, cfg.n_runs + 1))

    start_ps     = args.start     if args.start     is not None else cfg.extract_start_ps
    frequency_ps = args.frequency if args.frequency is not None else cfg.extract_frequency_ps

    if not setup.check_required_files():
        return

    if not args.skip_topology_test:
        setup.test_topology()

    if args.end is not None:
        end_ps = args.end
    elif cfg.extract_end_ps:
        end_ps = cfg.extract_end_ps
    else:
        first_tpr = Path(cfg.output_dir) / f"EQ_{run_ids[0]}" / "prod" / "prod.tpr"
        end_ps = get_production_length_ps(str(first_tpr), cfg.gmx)

    time_points_ns = frame_times_ns(start_ps, frequency_ps, end_ps)
    logger.info(
        f"Extracting {len(time_points_ns)} frames per run: "
        f"{time_points_ns[0]}–{time_points_ns[-1]} ns"
    )

    setup.create_inputs(run_ids, time_points_ns)
    logger.info("create_NE_NP_files.py finished.")


if __name__ == "__main__":
    main()
