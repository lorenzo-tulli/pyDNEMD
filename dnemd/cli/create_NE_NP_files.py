#!/usr/bin/env python3

"""
1. Build and test the perturbed topology (protein + solvent, no ligand).
2. Extract frames from each equilibrium run at the requested interval.
3. Create GROMACS input files for Non-Equilibrium (NE) simulations
   (perturbed topology, velocities NOT reassigned).
4. Create GROMACS input files for Null-Perturbation (NP) simulations
   (perturbed topology, velocities reassigned from Maxwell-Boltzmann distribution).

Directory layout created:
    output/
    ├── perturbed_topology/
    │   ├── perturbed_em.gro
    │   ├── perturbed_index.ndx
    │   └── topolperturb.top        ← user must provide
    ├── TRJDUMP_NE/
    │   └── NE_<run>/<time>ns/
    │       ├── <time>ns_NE.gro
    │       ├── Prod_RunNE.mdp
    │       └── MD_NE.tpr
    └── TRJDUMP_NP/
        └── NP_<run>/<time>ns/
            ├── <time>ns_NP.gro
            ├── Prod_RunNP.mdp
            └── MD_NP.tpr

The user must provide the file topolperturb.top by removing ligand entries from topol.top.

    topol.top:
    ...
    [ molecules ]
    ; Compound        #mols
    Protein             1
    Ligand              1
    SOL             39553
    NA                119
    CL                114

    topolperturb.top:
    ...
    [ molecules ]
    ; Compound        #mols
    Protein             1
    SOL             39553
    NA                119
    CL                114


Usage:
    python scripts/create_NE_NP_files.py --config examples/config.yaml \\
        --start 50000 --frequency 5000

    --start      : first frame to extract, in ps  (e.g. 50000 = 50 ns)
    --frequency  : interval between frames, in ps (e.g. 5000  =  5 ns)
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
    parser.add_argument("--start",     type=int,   required=True)
    parser.add_argument("--frequency", type=int,   required=True)
    parser.add_argument("--end",       type=int,   default=None)
    parser.add_argument("--run",       type=int,   default=None)
    parser.add_argument("--skip-topology-build", action="store_true")
    parser.add_argument("--skip-topology-test",  action="store_true")
    args = parser.parse_args()

    cfg     = Config.from_yaml(args.config)
    setup   = NESetup(cfg)
    run_ids = [args.run] if args.run else list(range(1, cfg.n_runs + 1))

    if not args.skip_topology_build:
        setup.build_topology()
    if not args.skip_topology_test:
        setup.test_topology()

    topolperturb = setup.perturb_dir / "topolperturb.top"
    if not topolperturb.exists():
        logger.error(
            f"topolperturb.top not found in {setup.perturb_dir}.\n"
            "Create it, then re-run with --skip-topology-build --skip-topology-test."
        )
        return

    if args.end is None:
        first_tpr = Path(cfg.output_dir) / f"EQ_{run_ids[0]}" / "prod" / "prod.tpr"
        end_ps = get_production_length_ps(str(first_tpr), cfg.gmx)
    else:
        end_ps = args.end

    time_points_ns = frame_times_ns(args.start, args.frequency, end_ps)
    logger.info(
        f"Extracting {len(time_points_ns)} frames per run: "
        f"{time_points_ns[0]}–{time_points_ns[-1]} ns"
    )

    setup.create_inputs(run_ids, time_points_ns)
    logger.info("create_NE_NP_files.py finished.")


if __name__ == "__main__":
    main()
